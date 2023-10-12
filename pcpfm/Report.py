from fpdf import FPDF
import datetime
import textwrap
import os
import subprocess
import platform

HEADER = 'PCPFM Report - '

class ReportPDF(FPDF):
    def header(self):
        global HEADER
        self.set_font('Arial', 'B', 15)
        self.cell(80)
        self.cell(30, 10, HEADER, 0, 0, 'C')
        self.ln(10)

    def footer(self):
        self.set_y(15)
        self.set_font('Arial', 'I', 8)

class Report():
    def __init__(self, experiment, parameters) -> None:
        self.experiment = experiment
        self.parameters = parameters

        global HEADER
        HEADER += self.experiment.experiment_directory.split("/")[-1]
        self.default_font = ['Arial', '', 12]
        self.report = self.initialize_report()
        self.max_width = round(self.report.line_width * 1000,0)
        self.style = self.preprocess_style(self.paramaters["requirements_txt"])

    def initialize_report(self):
        report = ReportPDF()
        report.add_page()
        report.set_font(self.default_font[0], self.default_font[1], self.default_font[2])
        return report
    
    @staticmethod
    def preprocess_style(style):
        for section in style:
            if "text" in section:
                section["text"] = style["texts"][section["text"]]
            else:
                section["text"] = None
        return section

    def create_report(self):
        for section in self.style:
            self.__getattribute__(section["section"], section)

    def all_TICs(self, section_desc):
        for acquisition in self.experiment.acquisitions:
            try:
                tic_path = acquisition.TICz()
                self.report.image(tic_path, w=self.max_width)
            except:
                pass

    def reset_font(self):
        self.report.set_font(self.default_font[0], self.default_font[1], self.default_font[2])

    def section_head(self, title):
        self.report.cell(80)
        self.report.set_font(self.default_font[0], 'B', self.default_font[2])
        self.report.cell(30, 10, title, 0, 0, 'C')
        self.reset_font()
        self.report.ln(5)

    def subsection_head(self, title):
        self.report.cell(80)
        self.report.cell(30, 10, title, 0, 0, 'C')
        self.report.ln(5)

    def section_line(self, content, options=None):
        options = set(options) if options is not None else set()
        if "bold" in options:
            self.report.set_font(self.default_font[0], 'B', self.default_font[2])
            self.report.cell(30, 10, content, 0, 0, "B")
            self.reset_font()
        else:
            self.report.cell(30, 10, content, 0, 0)
        self.report.ln(5)

    def section_text(self, text, options=None):
        text = ' '.join(text.split(None))
        text = ' '.join(text.split("\n"))
        for txt in textwrap.wrap(text, self.max_width / 2):
            self.section_line(txt, options=options)
        self.report.ln(5)

    def end_section(self):
        self.report.ln(5)

    def experiment_summary(self, text=None):
        pass

    def annotation_summary(self, section_desc):
        self.section_head("Annotation Summary")
        if 'text' in section_desc: 
            self.section_text(section_desc['text'])
        self.subsection_head("Feature Tables")
        self.section_line("Table Name, # Features, # MS1 Annotated Features, # MS2 Annotated Features", options=["bold"])

        for table in self.experiment.feature_tables.keys():
            try:
                feature_table = self.experiment.retrieve(table, True, False, True)
            except:
                feature_table = None
            if feature_table:
                num_features = feature_table.num_features
                num_ms2_annotations = 0
                if "MS2_annotations" in feature_table.feature_table.columns:
                    for x in feature_table.feature_table["MS2_annotations"]:
                        try:
                            x = str(x)
                            if 'score' in x:
                                num_ms2_annotations += 1
                        except:
                            pass
                else:
                    num_ms2_annotations = 0

                if "MS1_annotations" in feature_table.feature_table.columns:
                    ms1_annotated_features = [x for x in feature_table.feature_table["MS1_annotations"] if x and x != '[]']
                    num_ms1_annotations = len(ms1_annotated_features)
                else:
                    num_ms1_annotations = 0

                self.section_line(", ".join([str(x) for x in [table, num_features, num_ms1_annotations, num_ms2_annotations]]))
        self.section_head("")
        self.subsection_head("Empirical Compounds")
        self.section_line("empCpd Name, # Khipus, # MS1 Annotated Khipus, # MS2 Annotated Khipus", options=["bold"])
        for empcpd in self.experiment.empCpds.keys():
            try:
                empcpd_object = self.experiment.retrieve(empcpd, False, True, True)
            except:
                empcpd_object = None
            if empcpd_object:
                num_annotated_ms1 = 0
                num_annotated_ms2 = 0
                total = 0
                for kp_id, khipu in empcpd_object.dict_empCpds.items():
                    is_MS2_annotated = False
                    if "mz_only_db_matches" in khipu and khipu["mz_only_db_matches"]:
                        num_annotated_ms1 += 1
                    if "MS2_Spectra" in khipu and khipu["MS2_Spectra"]:
                        for spectrum in khipu["MS2_Spectra"]:
                            if "Annotations" in spectrum:
                                is_MS2_annotated = True
                    if is_MS2_annotated:
                        num_annotated_ms2 += 1
                    total += 1
                self.section_line(", ".join([str(x) for x in [empcpd, total, num_annotated_ms1, num_annotated_ms2]]))

    def table_summary(self, section_desc):
        self.section_head("Feature Table Summary")
        if 'text' in section_desc: 
            self.section_text(section_desc['text'])
        self.section_line("Table Name, Num Samples, Num Features", options=["bold"])
        for table in self.experiment.feature_tables.keys():
            try:
                feature_table = self.experiment.retrieve(table, True, False, True)
                self.section_line(", ".join([str(x) for x in [table, feature_table.num_samples, feature_table.num_features]]))
            except:
                pass

    def empcpd_summary(self, section_desc):
        self.section_head("empCpd Table Summary")
        if 'text' in section_desc: 
            self.section_text(section_desc['text'])
        self.section_line("EmpCpd Name, Num Khipus, Num Features", options=["bold"])
        for empcpd in self.experiment.empCpds.keys():
            try:
                empcpd_object = self.experiment.retrieve(empcpd, False, True, True)
                self.section_line(", ".join([str(x) for x in [empcpd, empcpd_object.num_khipus, empcpd_object.num_features]]))
            except:
                pass
            
    def command_history(self, section_desc):
        self.section_head("Command History")
        for command in self.experiment.command_history:
            self.section_line(command)

    def version_summary(self, section_desc):
        self.section_head("Software Version Summary")
        with open(self.parameters["requirements_txt"]) as req:
            for line in req:
                line = line.strip()
                output = subprocess.run(["pip", "show", line.rstrip()], capture_output=True, text=True)
                version = [x for x in output.stdout.split("\n") if x.startswith("Version")][0].split(": ")[-1].strip()
                self.section_line(":".join([line, version]))
            try:
                output = subprocess.run(["pip", "show", "pcpfm"], capture_output=True, text=True)
                version = [x for x in output.stdout.split("\n") if x.startswith("Version")][0].split(": ")[-1].strip()
                self.section_line(": ".join(["pcpfm", version]))
            except:
                pass
        self.section_line()
        self.section_line("OS: ", platform.system())
        self.section_line("Python Version: " + platform.python_version())
        self.section_line("Architecture: " + platform.machine())
        self.section_line("Uname: " + " ".join(platform.uname()))


    def timestamp(self, section_desc):
        timestamp_string = 'Report generated on ' + str(datetime.datetime.now())
        self.section_head("Timestamp")
        self.section_line(timestamp_string)

    def save(self, section_desc):
        if not section_desc["report_name"].endswith(".pdf"):
            section_desc["report_name"] = section_desc["report_name"] + ".pdf"
        out_path = os.path.join(os.path.abspath(self.experiment.experiment_directory), "reports", section_desc["report_name"])
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        self.report.output(os.path.join(os.path.abspath(self.experiment.experiment_directory), "reports", section_desc["report_name"]))

    def figure(self, section_desc):
        self.report.add_page()
        self.section_line("Table: " + section_desc["table"])
        self.report.ln(10)
        figure_path = os.path.join(os.path.abspath(self.experiment.experiment_directory), "QAQC_figs/" + section_desc["table"], section_desc["name"])
        self.report.image(figure_path, w=self.max_width)