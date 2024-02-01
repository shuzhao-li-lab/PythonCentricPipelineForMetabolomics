import matchms
import os 
import scipy.stats
import numpy as np
import shutil
from functools import partial

def row_to_dict(row, columns):
    _d = {}
    for c in columns:
        c = c.strip() if type(c) is str else c
        v = row[c].strip() if type(row[c]) is str else row[c]
        _d[c] = v
    return _d

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
        return matchms.similarity.__getattribute__(method_name)
    except:
        raise Exception("no matching similarity method named: ", method_name)

def extract_MS2_spectra(ms2_files, skip_filters=False):
    MS2_spectral_registry = {}
    for ms2_file in ms2_files:
        found = 0
        for spectrum in get_parser(ms2_file.split(".")[-1])(ms2_file, metadata_harmonization=True):
            found += 1
            if spectrum:
                MS2_spec_obj = process_ms2_spectrum(spectrum, filename=ms2_file, skip_filters=skip_filters)
                if MS2_spec_obj:
                    MS2_spectral_registry[MS2_spec_obj.id] = MS2_spec_obj
    return MS2_spectral_registry

def lazy_extract_MS2_spectra(ms2_files, mz_tree=None):
    if type(ms2_files) is str:
        ms2_files = [ms2_files]
    for ms2_file in ms2_files:
        for spectrum in get_parser(ms2_file.split(".")[-1])(ms2_file, metadata_harmonization=True):
            spectrum = matchms.filtering.add_precursor_mz(spectrum)
            if spectrum:
                try:
                    precursor_mz = spectrum.get('precursor_mz')
                except:
                    precursor_mz = None
            if precursor_mz and mz_tree:
                if mz_tree.at(precursor_mz):
                    try:
                        MS2_spec_obj = process_ms2_spectrum(spectrum, filename=ms2_file)
                        if MS2_spec_obj:
                            yield MS2_spec_obj
                    except:
                        pass
            elif mz_tree is None:
                MS2_spec_obj = process_ms2_spectrum(spectrum, filename=ms2_file)
                if MS2_spec_obj:
                    yield MS2_spec_obj

def process_ms2_spectrum(spectrum, filename="not_specified", min_peaks=3, skip_meta=False, skip_filters=False):
    """
    This is the default MS2 processing used by the pipeline. 

    This adds the precursor_mz, applies the default filters, normalizes
    the intensiies requires a minimum number of peaks and then converts
    the retention time to seconds. 

    :param filename: used to record what file the spectrum came from
    :param min_peaks: the min number of peaks a spectrum must have

    :return: dictionary summarize MS2 spectrum.
    """

    from .MSnSpectrum import MS2Spectrum

    if not skip_meta:
        spectrum = matchms.filtering.add_precursor_mz(spectrum)
        if not skip_filters:
            spectrum = matchms.filtering.default_filters(spectrum)
            spectrum = matchms.filtering.normalize_intensities(spectrum)
            spectrum = matchms.filtering.require_minimum_number_of_peaks(spectrum, min_peaks)

        if spectrum is None:
            return None
        try:
            spectrum.set('retention_time', spectrum.metadata['scan_start_time'][0] * 60)
        except:
            pass
        
        id = []
        if spectrum.get('precursor_mz'):
            id.append(str(float(spectrum.get('precursor_mz', ''))))
        if spectrum.get('retention_time'):
            id.append(str(float(spectrum.get('retention_time', ''))))
        id.append(filename)
        id = '_'.join(id)    
        try:
            reference_id = spectrum.get("compound_name")
        except:
            reference_id = id
        spectrum = MS2Spectrum(
                    id = id,
                    precursor_mz=float(spectrum.get('precursor_mz')),
                    precursor_rt=float(spectrum.get('retention_time')) if spectrum.get('retention_time') else None,
                    list_mz=[],
                    list_intensity=[],
                    matchms_spectrum=spectrum,
                    source=filename,
                    instrument="Unknown",
                    collision_energy="Unknown",
                    compound_name=reference_id)
    else:
        id = []
        id.append(str(float(spectrum['prec_mz'])))
        id.append(str(float(spectrum['rt'])))
        id.append(filename)
        id = '_'.join(id)
        precursor_mz = spectrum['prec_mz']
        precursor_rt = spectrum['rt']
        reference_id = spectrum['cpd_name']
        spectrum = spectrum['spectrum']
        if not skip_filters:
            spectrum = matchms.filtering.default_filters(spectrum)
            spectrum = matchms.filtering.normalize_intensities(spectrum)
            spectrum = matchms.filtering.require_minimum_number_of_peaks(spectrum, min_peaks)
        spectrum.set('retention_time', precursor_rt)
        spectrum.set('precursor_mz', precursor_mz)

        spectrum = MS2Spectrum(
                            id = id,
                            precursor_mz=float(precursor_mz),
                            precursor_rt=float(precursor_rt),
                            list_mz=[],
                            list_intensity=[],
                            matchms_spectrum=spectrum,
                            source=filename,
                            instrument="Unknown",
                            collision_energy="Unknown",
                            compound_name=reference_id)

    return spectrum
    
def search_for_mzml(sdir):
    """
    Given a directory, recusrively search for all files with an .mzml
    extension and return their paths.

    :param sdir: the directory in which to search

    :return: list of absolute mzML filepaths in sdir
    """
    mzml_found = []
    for d, _, fs in os.walk(sdir):
        for f in [f for f in fs if f.endswith(".mzML")]:
            mzml_found.append(os.path.join(os.path.abspath(d), f))
    return mzml_found


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
        elif hasattr(to_encode, "serialize"):
            return recursive_encoder(to_encode.serialize())
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

