'''
Usage:
  main.py assemble_experiment_from_CSV <experiment_directory> <experiment_name> <sequence_csv> <metadata_csv> [--filter=<filter>]
  main.py convert_to_mzML_local <experiment_directory> <mono_path> <ThermoRawFileConverter.exe_path> [--multi]
  main.py spectral_QCQA <experiment_directory> <standards_csv> <mz_search_tolerance_ppm> <rt_search_tolerance> <null_cutoff_percentile> <min_intensity> [--multi]
  main.py asari <experiment_directory>
'''

from docopt import docopt
import os
import logging
from Experiment import Experiment
import multiprocessing as mp
import subprocess


def job(job_desc):
    acquisition, spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory = job_desc
    acquisition.generate_report(spikeins, mz_ppm, rt_error, null_perc, min_intensity, text_report, output_directory)

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
        experiment = Experiment.construct_experiment_from_CSV(args['<experiment_name>'], args['<experiment_directory>'], args['<sequence_csv>'], args['<metadata_csv>'], filter=args['--filter'])
        experiment.save_experiment()
    if args['convert_to_mzML_local']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        experiment.convert_raw_to_mzML(args['<mono_path>'], args['<ThermoRawFileConverter.exe_path>'], multiprocessing=args['--multi'])
        experiment.save_experiment()
    if args['spectral_QCQA']:
        experiment = Experiment.load_experiment(os.path.join(args['<experiment_directory>'], "experiment.pickle"))
        # todo - fix this to take external standards
        spikeins = [
        # name, M+H, RT,
        #in extraction buffer
        #('13C6-D-glucose', 187.0908, 20),
        ('trimethyl-13C3-caffeine', 198.0977, None),
        ('15N-13C5-methionine', 156.0721, None),
        ('13C5-L-glutamate', 153.0722, None),
        ('15N2-uracil', 115.0286, None),
        ('15N-L-tyrosine', 183.0782, None),
        
        #spiked into samples directly
        ('15:0-18:1(d7) PC', 753.6134, None),
        ('18:1(d7) Lyso PC', 529.3994, None),
        ('15:0-18:1(d7) PE', 711.5664, None),
        ('18:1(d7) Lyso PE', 487.3524, None),
        ('15:0-18:1(d7) PG', 759.5875, None),
        ('15:0-18:1(d7) PI', 847.6036, None),
        ('15:0-18:1(d7) PS', 755.5562, None),
        ('15:0-18:1(d7)-15:0 TAG', 829.7985, None),
        ('15:0-18:1(d7) DAG', 605.5844, None),
        ('18:1(d7) MAG [H+]', 364.3429, None),
        ('18:1(d7) MAG [NH4+]', 381.3704, None),
        ('18:1(d7) Chol Ester', 675.6779, None),
        ('d18:1-18:1(d9) SM', 738.647, None),
        ('Cholesterol-d7', 411.4326, None)   
        ]

        if args['--multi']:
            job_descs = [(
                acquisition, 
                spikeins, 
                float(args["<mz_search_tolerance_ppm>"]), 
                float(args["<rt_search_tolerance>"]), 
                int(args["<null_cutoff_percentile>"]), 
                int(args["<min_intensity>"]), 
                True, 
                os.path.join(experiment.experiment_directory, "reports")) for acquisition in experiment.acquisitions]

            pool = mp.Pool(mp.cpu_count())
            pool.map(job, job_descs)
        else:
            for acquisition in experiment.acquisitions:
                print(acquisition.name)
                acquisition.generate_report(spikeins, float(args["<mz_search_tolerance_ppm>"]), float(args["<rt_search_tolerance>"]), int(args["<null_cutoff_percentile>"]), int(args["<min_intensity>"]), text_report=True, output_directory=os.path.join(experiment.experiment_directory, "reports"))
        experiment.save_experiment()

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

    