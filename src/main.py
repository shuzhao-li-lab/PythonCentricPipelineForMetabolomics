'''
Usage:
  main.py assemble_experiment_from_CSV <experiment_directory> <experiment_name> <sequence_csv> <metadata_csv>
  main.py convert_to_mzML_local <experiment_directory> <mono_path> <ThermoRawFileConverter.exe_path> [--multi]
  main.py asari_preprocessing <experiment_directory>
  main.py  
  main.py process_from_raw <experiment_directory> <experiment_name> <csv_path> <mono_path> <ThermoRawFileConverter.exe_path> [--multi]
'''

from docopt import docopt
import os
import logging
from asari.analyze import ext_estimate_min_peak_height
from Experiment import Experiment

def main(args):
    if '<experiment_directory>' in args:
        if os.path.isdir(args['<experiment_directory>']):
            pass
        else:
            os.makedirs(args['<experiment_directory>'])
        logging.basicConfig(filename = args['<experiment_directory>'] + "./analysis.log",
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s [%(filename)s:%(lineno)d] %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
        logger = logging.getLogger(__name__)
    if args['assemble_experiment_from_CSV']:
        experiment = Experiment.construct_experiment_from_CSV(args['<experiment_name>'], args['<experiment_directory>'], args['<sequence_csv>'], args['<metadata_csv>'])
        experiment.save_experiment()
    if args['convert_to_mzML_local']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        experiment.convert_raw_to_mzML(args['<mono_path>'], args['<ThermoRawFileConverter.exe_path>'], multiprocessing=args['--multi'])
        experiment.save_experiment()
    if args['asari_preprocessing']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        asari_params = {
            'mode': 'pos',
            'mz_tolerance_ppm': 5,
            'project_name': experiment.experiment_name,
            'outdir': experiment.asari_output,
            'input': experiment.converted_acquisitions_directory
        }
        params = ext_estimate_min_peak_height(os.listdir(experiment.converted_acquisitions_directory))
        print(params)

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

    