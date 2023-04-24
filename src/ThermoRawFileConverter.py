import logging
import subprocess
import multiprocessing as mp
import os
from abc import ABC, abstractmethod

class AbstractThermoRawFileConverter(ABC):
    log = logging.getLogger(__name__)

    @abstractmethod
    def convert(self, raw_acquisition, output_directory):
        return NotImplemented
    
    @abstractmethod
    def convert_multi(self, raw_acquisitions, output_directory):
        return NotImplemented
    

class ThermoRawFileConverter(AbstractThermoRawFileConverter):
    log = logging.getLogger(__name__)

    def __init__(self, mono_path, exe_path):
        self.log.info("Creating ThermoRawFileConverter for local install")
        self.log.info("\tMono path = " + mono_path)
        self.mono_path = mono_path
        self.log.info("\tThermoRawFileConverter executable path " + exe_path)
        self.exe_path = exe_path 

    def convert(self, raw_acquisition, output_directory):
        try:
            self.log.info("Converting " + raw_acquisition.raw_filepath)
            output_filepath = os.path.join(output_directory, os.path.basename(raw_acquisition.raw_filepath)).replace(".raw", ".mzML")
            print(output_filepath)
            subprocess.run([self.mono_path, self.exe_path, '-f=1', '-i', raw_acquisition.raw_filepath, '-b', output_filepath])
            self.log.info("Succesfully converted " + raw_acquisition.raw_filepath)
        except:
            self.log.exception("Error converting " + raw_acquisition.raw_filepath)

    def convert_multi(self, raw_acquisitions, output_directory):
        workers = mp.Pool(max(mp.cpu_count()-1, 1))
        jobs = []
        for raw_acquisition in raw_acquisitions:
            output_filepath = os.path.join(output_directory, os.path.basename(raw_acquisition.raw_filepath)).replace(".raw", ".mzML")
            jobs.append([self.mono_path, self.exe_path, '-f=1', '-i', raw_acquisition.raw_filepath, '-b', output_filepath])
        converted_filepaths = []
        for completed_process, raw_acquisition in zip(workers.map(subprocess.run, jobs), raw_acquisitions):
            if completed_process.returncode == 0:
                converted_filepaths.append(completed_process.args[6])
            else:
                converted_filepaths.append(None)
        return converted_filepaths