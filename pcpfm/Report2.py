import os
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple, Any
from fpdf import FPDF
import matplotlib.pyplot as plt

# ---------- Theme ----------
@dataclass(frozen=True)
class ReportTheme:
    font_family: str = "Arial"
    font_style: str = ""
    font_size: int = 12
    dpi: int = 300
    palette: str = "Okabe_Ito"  # hook for downstream fig styling

    @property
    def font_tuple(self) -> Tuple[str, str, int]:
        return (self.font_family, self.font_style, self.font_size)


# ---------- Experiment Mixins --------- #

def retrieve_figure(experiment, feature_table):
    params = experiment.parameters 




# ---------- PDF Wrapper ----------
class ReportPDF(FPDF):
    def __init__(self, header_text: str, theme: ReportTheme):
        self.header_text = header_text
        self.theme = theme
        super().__init__()

    def header(self):
        self.set_font(self.theme.font_family, 'B', 15)
        self.cell(80)
        self.cell(30, 10, self.header_text, 0, 0, 'C')
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font(self.theme.font_family, 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')


# ---------- Figure Handling ----------
class FigureRenderer:
    def __init__(self, dpi: int):
        self.dpi = dpi

    def __call__(self, fig, path: str):
        fig.savefig(path, dpi=self.dpi, bbox_inches='tight')
        plt.close(fig)


# ---------- Section Dispatcher ----------
class SectionDispatcher:
    def __init__(self, pdf: ReportPDF, theme: ReportTheme):
        self.pdf = pdf
        self.theme = theme
        self.renderer = FigureRenderer(theme.dpi)
        self.dispatch: Dict[str, Callable[[Dict[str, Any]], None]] = {
            'text': self._text,
            'figure': self._figure,
            'toc': self._toc,
            'cover': self._cover,
        }
        self._toc_positions: List[Tuple[str, int]] = []

    def run(self, section: Dict[str, Any]):
        kind = section.get('section')
        fn = self.dispatch.get(kind)
        if fn:
            fn(section)

    # --- handlers ---
    def _text(self, sec: Dict[str, Any]):
        self.pdf.set_font(*self.theme.font_tuple)
        title = sec.get('title')
        if title:
            self._section_title(title)
        txt = sec.get('text', '')
        self.pdf.multi_cell(0, 8, txt)

    def _figure(self, sec: Dict[str, Any]):
        fig = sec.get('figure')
        save_path = sec.get('save_path')
        caption = sec.get('caption', '')
        if fig and save_path:
            self.renderer(fig, save_path)
            self.pdf.image(save_path, w=180)
            if caption:
                self.pdf.ln(2)
                self.pdf.set_font(self.theme.font_family, 'I', 9)
                self.pdf.multi_cell(0, 5, caption)

    def _toc(self, _: Dict[str, Any]):
        # Insert TOC at end: capture current page, then render
        cur_y = self.pdf.get_y()
        self.pdf.add_page()
        self.pdf.set_font(self.theme.font_family, 'B', 16)
        self.pdf.cell(0, 10, 'Table of Contents', ln=1)
        self.pdf.set_font(self.theme.font_family, '', 12)
        for title, page in self._toc_positions:
            self.pdf.cell(0, 8, f'{title} ....... {page}', ln=1)
        self.pdf.set_y(cur_y)

    def _cover(self, sec: Dict[str, Any]):
        self.pdf.add_page()
        self.pdf.set_font(self.theme.font_family, 'B', 28)
        self.pdf.ln(60)
        title = sec.get('title', 'Report')
        self.pdf.cell(0, 20, title, 0, 1, 'C')
        subtitle = sec.get('subtitle', '')
        if subtitle:
            self.pdf.set_font(self.theme.font_family, '', 16)
            self.pdf.cell(0, 10, subtitle, 0, 1, 'C')

    # --- helpers ---
    def _section_title(self, title: str):
        self.pdf.set_font(self.theme.font_family, 'B', 14)
        self._toc_positions.append((title, self.pdf.page_no()))
        self.pdf.ln(6)
        self.pdf.cell(0, 10, title, ln=1)
        self.pdf.set_font(*self.theme.font_tuple)


# ---------- Style Validation ----------
class ReportStyle:
    def __init__(self, cfg: Dict[str, Any]):
        raw = cfg.get('sections', [])
        self.sections: List[Dict[str, Any]] = [s for s in raw if isinstance(s, dict) and 'section' in s]


# ---------- Main Report ----------
class Report:
    def __init__(self, experiment, parameters: Dict[str, Any], output_path: Optional[str] = None, theme: Optional[ReportTheme] = None):
        self.experiment = experiment
        self.parameters = parameters
        self.theme = theme or ReportTheme()

        header = f"PCPFM Report - {os.path.basename(experiment.experiment_directory)}"
        header = ''
        self.pdf = ReportPDF(header, self.theme)
        #self.pdf.add_page()
        self.pdf.set_font(*self.theme.font_tuple)

        style = ReportStyle(parameters['report_config'])
        dispatcher = SectionDispatcher(self.pdf, self.theme)

        for sec in style.sections:
            dispatcher.run(sec)

        self.output_path = output_path or os.path.join(self.experiment.experiment_directory, 'report.pdf')
        self.pdf.output(self.output_path)


# ---------- Example ----------
if __name__ == '__main__':
    import json
    import datetime
    from pcpfm.Experiment import Experiment

    exp = Experiment('P_HILICneg', '/Users/mitchjo/07202025_GPR32/Analysis/pcpfm/P_HILICneg')
    params = {
            "report_config": {
                "sections":
                    [
                        {'section': 'cover', 'subtitle': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")},
                        {'section': 'figure', 'table': 'full', 'name': 'pca', 'text': 'PCA'},
                        {'section': 'toc'}
                    ]
                }
            }
    Report(exp, params, output_path="./testing.pdf")
