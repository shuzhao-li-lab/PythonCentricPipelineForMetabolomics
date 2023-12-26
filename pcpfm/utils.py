import matchms
import os 
import scipy.stats
import numpy as np
import shutil

def get_parser(file_extension):
    """
    This will return the correct parser for a given MS2 file from 
    matchms.

    :param file_extension: the extension of a file, e.g., .msp, .mzml

    :return: appropriate parser method from matchms
    """
    try:
        return matchms.importing.__getattribute__("load_from_" + file_extension.lower())
    except:
        raise Exception("no matching parser for file type: ", file_extension)

def get_similarity_method(method_name):
    """
    This will return the specified similarity method from matchms.

    :param method_name: the name of a similarity method in matchms

    :return: appropriate similarity method from matchms
    """

    try:
        return matchms.similarity.__getattribute__(method_name)()
    except:
        raise Exception("no matching similarity method named: ", method_name)

def process_ms2_spectrum(spectrum, filename="not_specified", min_peaks=3):
    """
    This is the default MS2 processing used by the pipeline. 

    This adds the precursor_mz, applies the default filters, normalizes
    the intensiies requires a minimum number of peaks and then converts
    the retention time to seconds. 

    :param filename: used to record what file the spectrum came from
    :param min_peaks: the min number of peaks a spectrum must have

    :return: dictionary summarize MS2 spectrum.
    """

    try:
        spectrum = matchms.filtering.add_precursor_mz(spectrum)
        spectrum = matchms.filtering.default_filters(spectrum)
        spectrum = matchms.filtering.normalize_intensities(spectrum)
        spectrum = matchms.filtering.require_minimum_number_of_peaks(spectrum, min_peaks)
        spectrum.set('retention_time', spectrum.metadata['scan_start_time'][0] * 60)

        return {
            "exp_spectrum": spectrum, 
            "precursor_mz": float(spectrum.get('precursor_mz')),
            "precursor_rt": float(spectrum.get('retention_time')),
            "origin": filename,
            "annotations": []
        }
    except:
        return None
    
def search_for_mzml(sdir):
    """
    Given a directory, recusrively search for all files with an .mzml
    extension and return their paths.

    :param sdir: the directory in which to search

    :return: list of absolute mzML filepaths in sdir
    """
    return [os.path.join(os.path.abspath(d), f) for f, _, d in os.walk(sdir) if f.endwith(".mzML")]

def recursive_encoder(to_encode):
    """
    This method takes a datastructure and makes it JSON-ready by 
    recursively calling itself on the datastructure at all depths.
    This is needed because some objects such as Acquisition have their own
    method called JSON_repr that will encode it correctly whereas the
    JSON library will not. 

    Thus this can encode any datastructure comprised of other objeccs
    provided that they too are either a primitive type (e.g., str, int)
    or implement their own JSON_repr method. 

    :param to_encode: datastructure for conversion to JSON

    :return: JSON-encodable representation of to_encode
    """
    if to_encode is None:
        return None
    else:
        to_encode_type = type(to_encode)
        if to_encode_type in [str, int, float]:
            return to_encode
        elif to_encode_type in [list, set]:
            return [recursive_encoder(x) for x in to_encode]
        elif to_encode_type in [dict]:
            return {k: recursive_encoder(v) for k,v in to_encode.items()}
        elif hasattr(to_encode, "JSON_repr"):
            return to_encode.JSON_repr
        else:
            pass

# the following define mappings between strings and functions used in
# many places in the pcpfm. New modes can be added trivially by adding
# them to these mappings. 

correlation_modes = {
    "spearman": scipy.stats.spearmanr,
    "kendall": scipy.stats.kendalltau,
    "pearson": np.corrcoef
}

descriptive_stat_modes = {
    "median": np.median,
    "max": np.max,
    "min": np.min,
    "mean": np.mean
}

log_modes = {
    "log2": np.log2,
    "log10": np.log10
}

file_operations = {
    "link": os.link,
    "copy": shutil.copy,
    "move": shutil.move,
    "delete": shutil.rmtree,
}
