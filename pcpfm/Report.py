from fpdf import FPDF
import os
import json
from mass2chem.formula import calculate_mass, PROTON
import matplotlib.pyplot as plt
from . import FeatureTable

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
        self.report = ReportPDF('PCPFM Report - ' + self.experiment.experiment_directory.split("/")[-1])
        self.report.add_page()
        self.report.set_font(*self.default_font)
        self.max_width = round(self.report.line_width * 1000,0)
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
        top_level_fields = style.keys()
        for section in style["sections"]:
            for top_level_field in top_level_fields:
                if top_level_field in section and section[top_level_field] in style[top_level_field]:
                    section[top_level_field] = style[top_level_field][section[top_level_field]]
        valid_sections = []
        for section in style["sections"]:
            try:
                self.__getattribute__(section["section"])
                valid_sections.append(section)
            except:
                print(section["section"] + "is not a valid section!\nValid sections include: ")
                for method in dir(Report):
                    if not method.startswith("__") and not method.startswith("_Report__"):
                        print("\t", method)
        section_names = [x["section"] for x in valid_sections]
        if "save" not in section_names:
            raise Exception("Save section not found in report! will abort")
        return valid_sections

    def __create_report(self):
        """
        This function iterates through valid sections and generates the
        report element each section specifies.
        """
        for section in self.style:
            try:
                method = self.__getattribute__(section["section"])
                if "table" in section and section["table"] + "_cleaned" in self.experiment.feature_tables:
                    section["table"] = section["table"] + "_cleaned"
                    method(section)
                else:
                    method(section)
            except:
                print("Unable to processes section: \n", section)



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

    def __sub_section_head(self, title):
        """
        This writes a sub-section header. Large font, ala section_head
        but not in bold.

        TODO: this seems to be redundant.

        :param title: the title of the sub_section
        """
        self.report.cell(80)
        self.report.cell(30, 10, title, 0, 0, 'C')
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
        chunk_size = int(self.max_width // 2)
        while (i * chunk_size) < len(text):
            i += 1
            line = text[chunk_size * (i - 1): min(chunk_size * i, len(text))]
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

    def annotation_summary(self, section_desc):
        """
        This method summarizes the annotations generated during the 
        analysis on a per-table, per-empcpd basis. This counts the
        number of level-1, level-2, and level-4 annotations. 

        Requires: None
        """
        self.__section_head("Annotation Summary")
        if 'text' in section_desc: 
            self.__section_text(section_desc['text'])
        self.__sub_section_head("Feature Tables")
        self.__section_line("Table Name, # Features, # MS1 Annotated Features, # MS2 Annotated Features", options=["bold"])

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

                self.__section_line(", ".join([str(x) for x in [table, num_features, num_ms1_annotations, num_ms2_annotations]]))
        self.__section_head("")
        self.__sub_section_head("Empirical Compounds")
        self.__section_line("empCpd Name, # Khipus, # MS1 Annotated Khipus, # MS2 Annotated Khipus", options=["bold"])
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
                            if "Annotations" in spectrum and spectrum["Annotations"]:
                                is_MS2_annotated = True
                    if is_MS2_annotated:
                        num_annotated_ms2 += 1
                    total += 1
                self.__section_line(", ".join([str(x) for x in [empcpd, total, num_annotated_ms1, num_annotated_ms2]]))

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
            try:
                feature_table = self.experiment.retrieve(table, True, False, True)
                self.__section_line(", ".join([str(x) for x in [table, feature_table.num_samples, feature_table.num_features]]))
            except:
                pass

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
            try:
                empcpd_object = self.experiment.retrieve(empcpd, False, True, True)
                self.__section_line(", ".join([str(x) for x in [empcpd, empcpd_object.num_khipus, empcpd_object.num_features]]))
            except:
                pass
            
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

        import platform
        import subprocess
        self.__section_head("Software Version Summary")
        with open(self.parameters["requirements_txt"]) as req:
            for line in req:
                line = line.strip()
                output = subprocess.run(["pip", "show", line.rstrip()], capture_output=True, text=True)
                version = [x for x in output.stdout.split("\n") if x.startswith("Version")][0].split(": ")[-1].strip()
                self.__section_line(":".join([line, version]))
            try:
                output = subprocess.run(["pip", "show", "pcpfm"], capture_output=True, text=True)
                version = [x for x in output.stdout.split("\n") if x.startswith("Version")][0].split(": ")[-1].strip()
                self.__section_line(": ".join(["pcpfm", version]))
            except:
                pass
        self.__section_text("OS: ", platform.system())
        self.__section_text("Python Version: " + platform.python_version())
        self.__section_text("Architecture: " + platform.machine())
        self.__section_text("Uname: " + " ".join(platform.uname()))

    def computational_performance(self, section_desc):
        """
        This summarizes each command and computes the time required
        for that step. This is useful for benchmarking.

        Requires: None
        """
        import uuid
        import matplotlib.pyplot as plt
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
                pass
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
        import datetime
    
        timestamp_string = 'Report generated on ' + str(datetime.datetime.now())
        self.__section_head("Timestamp")
        self.__section_line(timestamp_string)

    def save(self, section_desc):
        """
        This saves the report pdf to the specified location

        Requires: "report_name"
        """
        if not section_desc["report_name"].endswith(".pdf"):
            section_desc["report_name"] = section_desc["report_name"] + ".pdf"
        out_path = os.path.join(os.path.abspath(self.experiment.experiment_directory), "reports", section_desc["report_name"])
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        self.report.output(os.path.join(os.path.abspath(self.experiment.experiment_directory), "reports", section_desc["report_name"]))

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
        
        figure_path = self.experiment.qaqc_figs + "/" + section_desc["table"] + "/_" + section_desc["name"] + ".png"      
        if os.path.exists(figure_path):
            self.report.add_page()
            self.__section_line("Table: " + section_desc["table"] + "  " + "Figure: " + section_desc["name"])
            self.report.ln(10)            
            self.report.image(figure_path, w=self.max_width)
        else:
            feature_table = self.experiment.retrieve(section_desc["table"], True, False, True)
            params_for_figure = {k: v for k,v in self.parameters.items()}
            params_for_figure['all'] = False
            params_for_figure['save_plots'] = True
            #params_for_figure['interactive_plots'] = False
            if section_desc["name"] in feature_table.qaqc_result_to_key:
                params_for_figure[feature_table.qaqc_result_to_key[section_desc["name"]]] = True
                feature_table.QAQC(params_for_figure)
                figure_path = self.experiment.qaqc_figs + "/" + section_desc["table"] + "/_" + section_desc["name"] + ".png"                
                self.report.add_page()
                self.__section_line("Table: " + section_desc["table"] + "  " + "Figure: " + section_desc["name"])
                self.report.ln(10)
                self.report.image(figure_path, w=self.max_width)

# this updates the docstring for 
getattr(Report, "figure").__doc__ += "\nValid name values are: \n\t" + "\n\t".join([x for x in FeatureTable.FeatureTable.qaqc_result_to_key.keys()]) 