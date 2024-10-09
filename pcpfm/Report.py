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
    This class is simply a wrapper around FPDF from fpdf library so
    that we can have a consistent header and footer for all pages

    :param header_text: text to be included on the header of each page
    """
    def __init__(self, header_text):
        self.header_text = header_text
        super().__init__()

    def header(self):
        """
        This function is called on every page of the report. It 
        generates the header for the page. 
        """
        self.set_font('Arial', 'B', 15)
        self.cell(80)
        self.cell(30, 10, self.header_text, 0, 0, 'C')
        self.ln(10)

    def footer(self):
        """
        This function is called on every page of the report. It 
        generates the footer. Currently this does nothing but it 
        is required. 
        """
        self.set_y(15)
        self.set_font('Arial', 'I', 8)

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
        self.report.add_page()
        self.report.set_font(*self.default_font)
        self.max_width = round(self.report.line_width * 900,0)
        self.style = self.__preprocess_style(self.parameters["report_config"])
        self.__create_report()

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
            try:
                if method:
                    if "table" in section:
                        if section["table"] + "_cleaned" in self.experiment.feature_tables:
                            section["table"] = section["table"] + "_cleaned"
                            method(section)
                        else:
                            method(section)
                    else:
                        method(section)
            except:
                pass
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

    def experiment_summary(self, section_desc):
        """
        This will list the empcpds and feature tables in the experiment.

        Requires: None
        """ 
        self.__section_head("Experiment Summary")
        if 'text' in section_desc:
            self.__section_text(section_desc['text'])
        self.__section_line("empCpd list", options=["bold"])
        for empcpd_moniker in self.experiment.empCpds.keys():
            self.__section_line(empcpd_moniker)
        self.__section_line(' ')
        self.__section_line("Feature Table list", options=["bold"])
        for feature_table_moniker in self.experiment.feature_tables.keys():
            self.__section_line(feature_table_moniker)


    def annotation_summary(self, section_desc):
        """
        This summarizes the annotations for each empcpd in the experiment.

        This counts the number of annotations per empcpd at each level.

        Args:
            Requires: None
        """     
        self.__section_head("Annotation Summary")
        if 'text' in section_desc:
            self.__section_text(section_desc['text'])
        self.__section_line("Name, #EmpCpds, #l4 Annotated, #l2 Annotated, #l1b Annotated, #l1a Annotated", options=["bold"])
        for empcpd_moniker in self.experiment.empCpds.keys():
            empcpds = self.experiment.retrieve_empCpds(empcpd_moniker, True)
            num_l4_annotated = 0
            num_l2_annotated = 0
            num_l1b_annotated = 0
            num_l1a_annotated = 0
            for kp in empcpds.dict_empcpds.values():
                num_l1b_annotated += int(bool(kp.get("Level_1b", [])))
                num_l4_annotated += int(bool(kp.get("Level_4")))
                has_l2 = False
                has_l1a = False
                for ms2_spectrum in kp.get("MS2_Spectra", []):
                    for annotation in ms2_spectrum.get("annotations", []):
                        if annotation["annotation_level"] == "Level_2":
                            has_l2 = True
                        elif annotation["annotation_level"] == "Level_1a":
                            has_l1a = True
            num_l2_annotated += int(has_l2)
            num_l1a_annotated += int(has_l1a)
            line = ", ".join([str(x) for x in [empcpd_moniker, 
                                               str(len(empcpds.dict_empcpds)), 
                                               num_l4_annotated, 
                                               num_l2_annotated, 
                                               num_l1b_annotated, 
                                               num_l1a_annotated]])
            self.__section_line(line)

    def table_summary(self, section_desc):
        """
        This summarizes the feature tables in the experiment. 

        Requires: None
        """

        self.__section_head("Feature Table Summary")
        if 'text' in section_desc:
            self.__section_text(section_desc['text'])
        self.__section_line("Table Name, Num Samples, Num Features", options=["bold"])
        tables = [x for x in self.experiment.feature_tables.keys() if "cleaned" not in x]
        for table in tables:
            ft = self.experiment.retrieve_feature_table(table, True)
            line = ", ".join([str(x) for x in [table, ft.num_samples, ft.num_features]])
            self.__section_line(line)

    def empcpd_summary(self, section_desc):
        """
        This summarizes the feature tables in the experiment. 

        Requires: None
        """
        self.__section_head("empCpd Table Summary")
        if 'text' in section_desc:
            self.__section_text(section_desc['text'])
        self.__section_line("EmpCpd Name, Num Khipus, Num Features", options=["bold"])
        for empcpd in self.experiment.empCpds.keys():
            empcpd_object = self.experiment.retrieve_empCpds(empcpd, True)
            self.__section_line(", ".join([str(x) for x in [empcpd, empcpd_object.num_khipus, empcpd_object.num_features]]))

    def command_history(self, section_desc):
        """
        This summarizes each command that has been executed in the 
        analysis. 

        Requires: None
        """

        self.__section_head("Command History")
        for command in self.experiment.command_history:
            self.__section_text(command)

    def version_summary(self, section_desc):
        """
        This summarizes each command that has been executed in the 
        analysis. 

        Requires: None
        """
        self.__section_head("Software Version Summary")
        self.__section_line(":".join(['pcpfm', pkg_resources.get_distribution('pcpfm').version]))
        _package = pkg_resources.working_set.by_key['pcpfm']
        for req in _package.requires():
            version = pkg_resources.get_distribution(req.name).version
            self.__section_line(":".join([req.name, version]))
        self.__section_line('')
        self.__section_line("OS: " + platform.system())
        self.__section_line("Python Version: " + platform.python_version())
        self.__section_line("Architecture: " + platform.machine())
        self.__section_text("Uname: " + " ".join(platform.uname()))

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
        timestamp_string = 'Report generated on ' + str(datetime.datetime.now())
        self.__section_head("Timestamp")
        self.__section_line(timestamp_string)

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
