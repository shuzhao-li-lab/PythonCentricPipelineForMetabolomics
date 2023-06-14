import subprocess
import multiprocessing as mp
import os

class ThermoRawFileConverter():
    def __init__(self, mono_path, exe_path):
        """
        This class wraps the ThermoRawFileParser for use with the pipeline.

        This really needs to be refactored or removed in the future.

        Args:
            mono_path (str): path to mono executable 
            exe_path (str): path to ThermoRawFileParser.exe
        """        
        self.mono_path = mono_path
        self.exe_path = exe_path 

    def convert(self, raw_acquisition, output_directory):
        """
        Convert a raw acquisition and put the resulting .mzML into the output_directory

        This is not used currently, not sure if this fully works. 

        Args:
            raw_acquisition (Acquisition): Acquisition object to convert
            output_directory (str): output directory
        """        
        try:
            output_filepath = os.path.join(output_directory, os.path.basename(raw_acquisition.raw_filepath)).replace(".raw", ".mzML")
            subprocess.run([self.mono_path, self.exe_path, '-f=1', '-i', raw_acquisition.raw_filepath, '-b', output_filepath, ' >/dev/null'])
        except:
            pass

    def convert_multi(self, raw_acquisitions, output_directory):
        """
        Convert multiple acquisitions and put the resulting .mzML into the output directory in parallel

        Args:
            raw_acquisitions (list[Acquisition]): list of Acquisition objects to convert
            output_directory (str): path to output directory

        Returns:
            list[str]: list of resulting .mzML filepaths co-indexed with the acquisitions.
        """        
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