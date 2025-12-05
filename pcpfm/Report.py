"""
This module implements the report object which is a wrapper around an fpdf object. 

The reports are meant to be a high-level overview of an experiment, the feature tables,
empcpds, etc and summarize some of the qaqc results.
"""

import os
import platform
import sys
import uuid
import datetime
import textwrap
from pip._vendor import pkg_resources
from fpdf import FPDF
import matplotlib.pyplot as plt
from . import FeatureTable
from . import utils

class ReportPDF(FPDF):
    """
    Wrapper around FPDF with consistent header and footer.
    """
    def __init__(self, header_text):
        self.header_text = header_text
        super().__init__()

    def header(self):
        print(self.page_no())
        if self.page_no() > 1:
            self.set_font('Arial', 'B', 15)
            #self.cell(0, 10, self.header_text, ln=True, align='C')
            self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        if self.page_no() > 1:
            self.cell(0, 10, f'Page {self.page_no() - 1} - {self.header_text}', 0, 0, 'C')

class Report():
    """
    The report object allows the creation of a pdf report. The reports 
    are defined using templates. 

    A template is a JSON-formatted list of dictionaries. Each dictionary 
    corresponds one-to-one with a method of this object. Which method
    the dictionary is intended to call is defined by the dictionary's
    "section" field. The section field value should be the name of a 
    method for this object. 
    """
    def __init__(self, experiment, parameters) -> None:
        self.experiment = experiment
        self.parameters = parameters

        self.default_font = ['Arial', '', 12]
        report_title = 'PCPFM Report - ' + self.experiment.experiment_directory.split("/")[-1]
        self.report = ReportPDF(report_title)
        #self.report.add_page()
        self.report.set_font(*self.default_font)
        self.max_width = round(self.report.line_width * 900,0)
        self.style = self.__preprocess_style(self.parameters["report_config"])
        self.__create_report()

    def table_of_contents(self, section_desc):
        """
        Adds a Table of Contents page. Assumes each section = 1 page.
        Automatically lists all sections before 'save', skipping 'cover_page' and TOC itself.
        """
        self.report.add_page()
        self.__section_head("Table of Contents")
        
        page_num = 2  # cover = 1, TOC = 2
        toc_lines = []
        for section in self.__preprocess_style(self.parameters["report_config"]):
            name = section["section"]
            if name in {"cover_page", "table_of_contents", "save"}:
                continue
            if section['section'] == 'figure':
                label = section['name']
            else:
                label = name.replace("_", " ").title()
            if "table" in section and section["table"] != "*":
                label += f" ({section['table']})"
            if (label, page_num) not in toc_lines:
                toc_lines.append((label, page_num))
                page_num += 1

        # Set uniform layout
        for label, page in toc_lines:
            label_str = f"{label}"
            label_width = self.report.get_string_width(label_str)
            dots = "." * int((180 - label_width - 10) / self.report.get_string_width("."))
            self.report.cell(0, 10, f"{label_str} {dots} {page}", ln=True)

    def __section_table(self, rows, col_widths=None, header=True, font_size=12):
        """
        Draws a well-aligned table in the PDF.
        
        Args:
            rows (list of lists): Each sublist is a row of strings/numbers.
            col_widths (list of floats): Optional fixed column widths in mm.
            header (bool): Whether the first row is a header.
            font_size (int): Font size for table content.
        """
        self.report.set_font(self.default_font[0], 'B' if header else '', font_size)
        
        # Auto-compute column widths if not provided
        if not col_widths:
            max_cols = max(len(r) for r in rows)
            total_width = self.report.w - self.report.l_margin - self.report.r_margin
            col_widths = [total_width / max_cols] * max_cols

        for i, row in enumerate(rows):
            if i == 1 and header:
                self.report.set_font(self.default_font[0], '', font_size)
            for j, cell in enumerate(row):
                self.report.cell(col_widths[j], 10, str(cell), border=1)
            self.report.ln()

    def __preprocess_style(self, style):
        """
        This function takes the provided style and checks that it has 
        the required save section and that the defined sections are 
        valid. Invalid sections are removed here. Furthermore, fields
        such as "texts" that can occur at the same level as "sections" 
        are then inserted into the approprate location. 

        :param style: the JSON configuration of the report

        :return: the style with text added and invalid methods deleted
        """
        for section in style["sections"]:
            for top_level_field in style.keys():
                if top_level_field in section:
                    if section[top_level_field] in style[top_level_field]:
                        section[top_level_field] = style[top_level_field][section[top_level_field]]
        valid_sections = []
        for section in style["sections"]:
            try:
                getattr(self, section["section"])
                valid_sections.append(section)
            except AttributeError:
                print(section["section"] + "is not a valid section!\nValid sections include: ")
                for method in dir(Report):
                    if not method.startswith("__") and not method.startswith("_Report__"):
                        print("\t", method)
        section_names = [x["section"] for x in valid_sections]
        if "save" not in section_names:
            print("Save section not found in report! will abort")
            sys.exit()

        expanded_sections_ft = []
        all_feature_tables = list(self.experiment.feature_tables.keys())
        for section in valid_sections:
            table_moniker = section.get("table", None)
            if table_moniker == "*":
                for feature_table in all_feature_tables:
                    new_section = dict(section)
                    new_section["table"] = feature_table
                    expanded_sections_ft.append(new_section)
            elif isinstance(table_moniker, list):
                for feature_table in table_moniker:
                    new_section = dict(section)
                    new_section["table"] = feature_table
                    expanded_sections_ft.append(new_section)
            else:
                expanded_sections_ft.append(dict(section))

        expanded_sections_ft_empcpd = []
        all_empcpds = list(self.experiment.empCpds.keys())
        for section in expanded_sections_ft:
            empcpd_moniker = section.get("empcpd", None)
            if empcpd_moniker == "*":
                for empcpd in all_empcpds:
                    new_section = dict(section)
                    new_section["empcpd"] = empcpd
                    expanded_sections_ft_empcpd.append(new_section)
            else:
                expanded_sections_ft_empcpd.append(dict(section))
        return expanded_sections_ft_empcpd

    def __create_report(self):
        """
        This function iterates through valid sections and generates the
        report element each section specifies.

        The way this works is that the section field should have a value whose name is the same
        as a method of the object, it finds and executes that method. As a result, every method
        must have the same call signature which results in some undesirable warnings from the IDE
        but is otherwise fine. This can be fixed in the future using something like the inspect
        module to see what the function requires and only pass that data.
        """
        for section in self.style:
            method = getattr(self, section["section"])

            if method:
                if "table" in section:
                    if section["table"] + "_cleaned" in self.experiment.feature_tables:
                        section["table"] = section["table"] + "_cleaned"
                        method(section)
                    else:
                        method(section)
                else:
                    method(section)

        self.experiment.save()


    def __reset_font(self):
        """
        This method resets the font of the report to the default. This 
        is needed because the report's font is set at a report-wide 
        level and must be changed each time a different font, size, or
        typesetting is required. 
        """
        self.report.set_font(*self.default_font)

    def __section_head(self, title):
        """
        This writes the header for a section. This uses a larger font
        and sets it to be in bold.
        
        :param title: the title of the section
        """
        self.report.cell(80)
        self.report.set_font(self.default_font[0], 'B', self.default_font[2])
        self.report.cell(30, 10, title, 0, 0, 'C')
        self.__reset_font()
        self.report.ln(5)

    def __section_line(self, content, options=None):
        """
        This writes a line in a section. Optionally, the font can be 
        made bold by passing 'B' in options. 

        TODO: all of these functions for writing lines can be 
        unified with a single command that takes multiple params

        :param content: what to write in the line
        :param options: an iterable with text options valid options
            include: 'B' for bold.
        """
        options = set(options) if options is not None else set()
        if "bold" in options:
            self.report.set_font(self.default_font[0], 'B', self.default_font[2])
            self.report.cell(30, 10, content, 0, 0, "B")
            self.__reset_font()
        else:
            self.report.cell(30, 10, content, 0, 0)
        self.report.ln(5)

    def __section_text(self, text, options=None):
        """
        This writes a block of text in a section, one line at a time
        using the __section_line function. This also handles chopping 
        up a line into pieces that can fit on the page (i.e., word 
        wrapping)

        TODO: all of these functions for writing lines can be 
        unified with a single command that takes multiple params

        :param text: what to write in the line
        :param options: an iterable with text options valid options
            include: 'B' for bold.
        """
        text = ' '.join(text.split(None))
        text = ' '.join(text.split("\n"))
        i = 0
        for line in textwrap.wrap(text, width=95):
            self.__section_line(line, options=options)
        self.report.ln(5)
        self.report.ln(5)

    def TICs(self, section_desc):
        """
        This generates, if not pre-existing, and includes the TIC of 
        each acquisition in the experiment to the report. 

        Requires: None
        """
        for acquisition in self.experiment.acquisitions:
            try:
                tic_path = acquisition.TIC()
                self.report.image(tic_path, w=self.max_width)
            except:
                print(section_desc)

    def study_summary(self, section_desc):
        pass

    def experiment_summary(self, section_desc):
        """
        Summarize key aspects of the experiment concisely.
        """
        self.report.add_page()
        self.__section_head("Experiment Overview")
        
        exp = self.experiment

        total_samples = len(exp.acquisitions)
        
        type_counts = {}
        for acq in exp.acquisitions:
            stype = getattr(acq, "sample_type", "Unknown")
            if stype == "Unknown":
                stype = "Study Sample"
            type_counts[stype] = type_counts.get(stype, 0) + 1
        type_summary = ", ".join([f"{count} {stype}" for stype, count in type_counts.items()])

        #num_metadata_tags = len(getattr(exp.metadata, "columns", {}))
        num_feature_tables = len(exp.feature_tables)

        summary_text = (
            f"A total of {total_samples} samples were analyzed.\n"
            f"Samples by type: {type_summary}.\n"
            #f"{num_metadata_tags} metadata tags were provided.\n"
            f"The experiment has {num_feature_tables} feature tables: ."
        )

        self.__section_text(summary_text)

    def annotation_summary(self, section_desc):
        """
        Summarizes annotations for each empCpd in a structured table.
        """
        self.report.add_page()
        self.__section_head("Annotation Summary")
        if 'text' in section_desc:
            self.__section_text(section_desc['text'])

        rows = [["Name", "#EmpCpds", "#L4 Annotated", "#L2 Annotated", "#L1b Annotated", "#L1a Annotated"]]

        for empcpd_moniker in self.experiment.empCpds.keys():
            try:
                empcpds = self.experiment.retrieve_empCpds(empcpd_moniker, True)
                num_l4, num_l2, num_l1b, num_l1a = 0, 0, 0, 0

                for kp in empcpds.dict_empcpds.values():
                    num_l1b += int(bool(kp.get("Level_1b", [])))
                    num_l4 += int(bool(kp.get("Level_4")))
                    annotations = [annotation["annotation_level"]
                                for spectrum in kp.get("MS2_Spectra", [])
                                for annotation in spectrum.get("annotations", [])]
                    num_l2 += int("Level_2" in annotations)
                    num_l1a += int("Level_1a" in annotations)

                rows.append([
                    empcpd_moniker,
                    len(empcpds.dict_empcpds),
                    num_l4,
                    num_l2,
                    num_l1b,
                    num_l1a
                ])
            except:
                pass

        self.__section_table(rows)


    def summary(self, section_desc):
        self.report.add_page()
        self.__section_head("Experiment Summary")
        subsections = [
            ("sample_desc_txt", "Sample Description"),
            ("experiment_goal_txt", "Goal"),
            ("assay_summary", "Assay Summary"), 
            ("analysis_summary", "Analysis Summary"),
            ("results_summary", "Results Summary")
        ]
        for (key, title) in subsections:
            self.__section_line(f"{title}")
            self.report.ln(5)
            self.__section_text(section_desc.get(key, ''))


    def table_summary(self, section_desc):
        """
        This summarizes the feature tables in the experiment. 

        Requires: None
        """
        self.report.add_page()
        self.__section_head("Feature Table Summary")
        if 'text' in section_desc:
            self.__section_text(section_desc['text'])
        tables = [x for x in self.experiment.feature_tables.keys() if "__internal" not in x]
        rows = [["Table Name", "Num Samples", "Num Features"]]
        for table in tables:
            ft = self.experiment.retrieve_feature_table(table, True)
            rows.append([table, ft.num_samples, ft.num_features])
        self.__section_table(rows)

    def empcpd_summary(self, section_desc):
        """
        Summarizes the empCpds in the experiment.
        """
        self.report.add_page()
        self.__section_head("empCpd Table Summary")
        if 'text' in section_desc:
            self.__section_text(section_desc['text'])

        rows = [["EmpCpd Name", "Num Khipus", "Num Features"]]
        for empcpd in self.experiment.empCpds.keys():
            try:
                obj = self.experiment.retrieve_empCpds(empcpd, True)
                rows.append([empcpd, obj.num_khipus, obj.num_features])
            except:
                pass
        
        self.__section_table(rows)

    def command_history(self, section_desc):
        """
        Summarizes each command that has been executed in the analysis.
        """
        self.report.add_page()
        self.__section_head("Command History")

        rows = [["Timestamp", "Command"]]
        for entry in self.experiment.command_history:
            ts, cmd = entry.split(":", 1)
            timestamp = datetime.datetime.fromtimestamp(float(ts)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append([timestamp, cmd])

        self.__section_table(rows, col_widths=[50, self.report.w - self.report.l_margin - self.report.r_margin - 50])

    def cover_page(self, section_desc):
        """
        Adds a cover page to the report.

        Optional fields in section_desc:
            - title (str): Main title of the report
            - subtitle (str): Subtitle or project description
            - author (str): Name of the author or lab
            - date (str): Override date string (default: today)
            - logo_path (str): Path to a logo image
        """
        self.report.add_page()
        self.report.set_font('Arial', 'B', 24)
        self.report.ln(100)
        if "title" in section_desc:
            self.report.cell(0, 40, section_desc["title"], ln=True, align='C')
        if "subtitle" in section_desc:
            self.report.set_font('Arial', '', 16)
            self.report.cell(0, 10, section_desc["subtitle"], ln=True, align='C')
        if "author" in section_desc:
            self.report.ln(10)
            self.report.set_font('Arial', '', 12)
            self.report.cell(0, 10, "Prepared by: " + section_desc["author"], ln=True, align='C')
        date_str = section_desc.get("date", datetime.datetime.today().strftime('%B %d, %Y    %H:%M:%S'))
        self.report.ln(5)
        self.report.cell(0, 10, "Date: " + date_str, ln=True, align='C')
        if "logo_path" in section_desc and os.path.exists(section_desc["logo_path"]):
            self.report.image(section_desc["logo_path"], x=80, w=50)

    def version_summary(self, section_desc):
        """
        This summarizes each command that has been executed in the 
        analysis. 

        Requires: None
        """
        self.report.add_page()
        self.__section_head("Software Metadata Summary")
        self.report.ln(5)
        __PCPFM_version = pkg_resources.get_distribution('pcpfm').version

        __version_text = f'PCPFM version {__PCPFM_version} was used for data processing\n'
        __version_text += 'The following software libraries were employed for the analysis: '

        submodules = []
        for req in pkg_resources.working_set.by_key['pcpfm'].requires():
            submodules.append(f'{req.name}v.{pkg_resources.get_distribution(req.name).version}')
        __version_text += " ,".join(submodules)
        self.__section_text(__version_text)
        self.report.ln(5)

        __machine_description = f'Report generated and analysis performed using a '
        __machine_description += f"{platform.machine()} system running {platform.system()} as the operating system."
        __uname_string = ",".join(platform.uname())
        __machine_description += f"Python Version was {platform.python_revision()} and system uname was {__uname_string}"
        self.__section_text(__machine_description)

    def computational_performance(self, section_desc):
        """
        This summarizes each command and computes the time required
        for that step. This is useful for benchmarking.

        Requires: None
        """

        self.__section_head("Computational Performance")
        tn_minus_one = None
        current_time = None
        command_order, time_required = [], []
        start_time = None
        for command in self.experiment.command_history:
            tn_minus_one = current_time
            current_time = float(command.split(":")[0])
            if command.endswith("start_analysis"):
                start_time = float(current_time)
            else:
                time_sums = (current_time - tn_minus_one) / 60
                command_order.append(command)
                time_required.append(time_sums)
        total_time = current_time - start_time
        self.__section_text("Total time for analysis: " + str(total_time / 60) + " minutes")
        plt.bar(command_order, time_required)
        plt.title("Time Required per Command")
        name = "/tmp/" + str(uuid.uuid4) + ".png"
        plt.savefig(name)
        self.report.image(name)

    def timestamp(self, section_desc):
        """
        This will timestamp the report.

        Requires: None
        """
        pass
        #timestamp_string = 'Report generated on ' + str(datetime.datetime.now())
        #self.__section_head("Timestamp")
        #self.__section_line(timestamp_string)

    def save(self, section_desc):
        """
        This saves the report pdf to the specified location

        Requires: "report_name"
        """
        output_subdir = os.path.abspath(self.experiment.output_subdirectory)
        report_path = os.path.join(output_subdir, section_desc["report_name"])
        if not report_path.endswith(".pdf"):
            report_path += ".pdf"
        print("saving: ", report_path)
        self.report.output(report_path)

    def figure_live(self, section_desc):
        feature_table = self.experiment.retrieve_feature_table(section_desc["table"], True)
        params_for_figure = dict(self.parameters)
        params_for_figure['all'] = False
        params_for_figure['save_plots'] = True
        feature_table.generate_figure_params(params_for_figure)
        figure_path = feature_table.save_fig_path(section_desc["name"])
        if section_desc["name"] in feature_table.qaqc_result_to_key:
                    params_for_figure[feature_table.qaqc_result_to_key[section_desc["name"]]] = True
                    feature_table.QAQC(params_for_figure)


    def figure(self, section_desc):
        """
        This inserts a figure into the report. Since all figures are 
        currently QAQC figures from feature tables, this method 
        requires specifying a table name and a qaqc result to 
        visualize. 

        Requires: "table" - the moniker for the table
                  "name" - name for the qaqc result to add figure of

        The valid fields of name are populated at runtime for this 
        method and will be displayed if an incorrect name field is 
        provided.

        """
        feature_table = self.experiment.retrieve_feature_table(section_desc["table"], True)
        params_for_figure = dict(self.parameters)
        params_for_figure['all'] = False
        params_for_figure['save_plots'] = True
        feature_table.generate_figure_params(params_for_figure)
        figure_path = feature_table.save_fig_path(section_desc["name"])
        try:
            if os.path.exists(figure_path):
                self.report.add_page()
                self.__section_line("Table: " + section_desc["table"] + "  " + "Figure: " + section_desc["name"])
                self.report.ln(10)
                self.report.image(figure_path, w=self.max_width-10)
            else:
                feature_table = self.experiment.retrieve_feature_table(section_desc["table"], True)
                params_for_figure = dict(self.parameters)
                params_for_figure['all'] = False
                params_for_figure['save_plots'] = True
                if section_desc["name"] in feature_table.qaqc_result_to_key:
                    params_for_figure[feature_table.qaqc_result_to_key[section_desc["name"]]] = True
                    feature_table.QAQC(params_for_figure)
                    self.report.add_page()
                    self.__section_line("Table: " + section_desc["table"] + "  " + "Figure: " + section_desc["name"])
                    self.report.ln(10)
                    self.report.image(figure_path, w=self.max_width)
            if "text" in section_desc and section_desc["text"]:
                self.__section_text(section_desc["text"])
        except:
            pass

# this updates the docstring for report
qaqc_names = list(FeatureTable.FeatureTable.qaqc_result_to_key.keys())
getattr(Report, "figure").__doc__ += "\nValid name values are: \n\t" + "\n\t".join(qaqc_names)
