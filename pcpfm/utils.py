"""
Misc helper functions.

"""
import os
import shutil
import sys
import matchms
import scipy.stats
import numpy as np
import pandas as pd

from mass2chem.formula import calculate_formula_mass
from matchms.Spectrum import Spectrum
from .MSnSpectrum import MS2Spectrum

def flatten_nested_dicts(nested_dicts):
    """
    For nested dictionaries, i.e., where the value for a key in a dictionary is a dictionary
    return a flat, dictionary where the nested keys become top level keys. 

    If the same key is used at different levels, the one encountered last will be the 
    top level key, value pair.

    Args:
        nested_dicts (dict): a dictionary that has a dictionary as a value

    Returns:
        dict: flattened dict
    """    
    _d = {}
    for k, v in nested_dicts.items():
        if not isinstance(v, dict):
            _d[k] = v
        else:
            _d.update(flatten_nested_dicts(v))
    return _d


def extract_CD_csv(
    cd_csvs, ionization_mode, min_peaks=1, precursor_to_use="Confirmed", lazy=True
):
    """
    For a list of compound discover (CD) CSV libraries, read them into a form that is 
    usable for level1a annotation.

    In lazy mode, this yields the library entries, else a list of them is return.

    Args:
        CD_csv (list): list of paths to cd csv libraries
        ionization_mode (str): the ionization mode for the experiment (and library)
        min_peaks (int, optional): entries with fewer than these many peaks are discarded. Defaults to 1.
        precursor_to_use (str, optional): CD library has "confirmed" and "extracted" precursors, specifies which to use as the precursor m/z. Defaults to "Confirmed".
        lazy (bool, optional): determines if this is a generator or returns a list. Defaults to True.

    Returns:
        list: list of MS2 spectra

    Yields:
        MSnSpectrum: object representing an MSnSpectrum
    """
    standards_spectra = []
    for cd_csv in cd_csvs:
        for standard in pd.read_csv(cd_csv).to_dict(orient="records"):
            mzs, intensities = [], []
            for i in range(10):
                mzs.append(
                    standard["Confirm Extracted." + str(i)]
                    if i and "Confirm Extracted." + str(i) in standard
                    else standard["Confirm Extracted"]
                )
                intensities.append(
                    standard["Target Ratio." + str(i)]
                    if i and "Target Ratio." + str(i) in standard
                    else standard["Target Ratio"]
                )
            mzs = [mz for mz in mzs if (mz and not np.isnan(mz))]
            intensities = [i for i in intensities if (i and not np.isnan(i))]
            mzs, intensities = zip(*sorted(zip(mzs, intensities)))
            if len(mzs) == len(intensities) and len(mzs) >= min_peaks:
                spectrum = Spectrum(
                    mz=np.asarray([float(x) for x in mzs]),
                    intensities=np.asarray([float(x) for x in intensities]),
                    metadata=standard,
                )
                # assume proton for ionization adduct
                if ionization_mode == "pos":
                    standard["Theoretical_Precursor"] = (
                        calculate_formula_mass(standard["ChemicalFormula"])
                        + calculate_formula_mass("H")
                        - calculate_formula_mass("e")
                    )
                else:
                    standard["Theoretical_Precursor"] = (
                        calculate_formula_mass(standard["ChemicalFormula"])
                        - calculate_formula_mass("H")
                        + calculate_formula_mass("e")
                    )
                spectrum = process_ms2_spectrum(
                    {
                        "rt": standard["RT"] * 60,
                        "prec_mz": standard["Confirm Precursor"]
                        if precursor_to_use == "Confirmed"
                        else standard["Theoretical_Precursor"],
                        "cpd_name": standard["CompoundName"],
                        "spectrum": spectrum,
                    },
                    filename=cd_csv,
                    min_peaks=min_peaks,
                    skip_meta=True,
                    skip_filters=False,
                )
            if spectrum and not lazy:
                standards_spectra.append(spectrum)
            if spectrum and lazy:
                yield spectrum
    return standards_spectra


def get_parser(file_extension):
    """
    This will return the correct parser for a given MS2 file from
    matchms.

    :param file_extension: the extension of a file, e.g., .msp, .mzml

    :return: appropriate parser method from matchms
    """
    method = getattr(matchms.importing, "load_from_" + file_extension.lower())
    if method:
        return method
    print("No matching parser found for: ", file_extension)
    sys.exit()


def get_similarity_method(method_name):
    """
    This will return the specified similarity method from matchms.

    :param method_name: the name of a similarity method in matchms

    :return: appropriate similarity method from matchms
    """

    method = getattr(matchms.similarity, method_name)
    if method:
        return method
    print("No matching similarity method found for: ", method_name)
    sys.exit()


def lazy_extract_ms2_spectra(ms2_files, mz_tree=None):
    """
    This method takes a list of ms2 files and yields the MS2 spectrum 
    for each entry. 

    Args:
        ms2_files (list): list of ms2_files, can be any format matchms supports.
        mz_tree (interval_tree, optional): interval tree of precursor mzs, if provided, only 
        spectra whose precuror matches the tree will be returned. Defaults to None.

    Yields:
        MSnSpectrum: object representing an MSnSpectrum
    """
    for ms2_file in [ms2_files] if isinstance(ms2_files, str) else ms2_files:
        for spectrum in get_parser(ms2_file.split(".")[-1])(
            ms2_file, metadata_harmonization=True
        ):
            spectrum = matchms.filtering.add_precursor_mz(spectrum)
            if spectrum:
                precursor_mz = spectrum.get("precursor_mz", None)
            if precursor_mz:
                if mz_tree:
                    if mz_tree.at(precursor_mz):
                        ms2_spec_obj = process_ms2_spectrum(spectrum, filename=ms2_file)
                        if ms2_spec_obj:
                            yield ms2_spec_obj
                else:
                    ms2_spec_obj = process_ms2_spectrum(spectrum, filename=ms2_file)
                    if ms2_spec_obj:
                        yield ms2_spec_obj


def process_ms2_spectrum(
    spectrum, filename="not_specified", min_peaks=1, skip_meta=False, skip_filters=False
):
    """
    This is the default MS2 processing used by the pipeline.

    This adds the precursor_mz, applies the default filters, normalizes
    the intensiies requires a minimum number of peaks and then converts
    the retention time to seconds.

    :param filename: used to record what file the spectrum came from
    :param min_peaks: the min number of peaks a spectrum must have

    :return: dictionary summarize MS2 spectrum.
    """

    if not skip_meta:
        spectrum = matchms.filtering.add_precursor_mz(spectrum)
        if not skip_filters:
            spectrum = matchms.filtering.default_filters(spectrum)
            spectrum = matchms.filtering.normalize_intensities(spectrum)
            spectrum = matchms.filtering.require_minimum_number_of_peaks(
                spectrum, min_peaks
            )

        if spectrum is None:
            return None
        try:
            spectrum.set("retention_time", spectrum.metadata["scan_start_time"][0] * 60)
        except:
            pass

        spec_id = []
        if spectrum.get("precursor_mz"):
            spec_id.append(str(float(spectrum.get("precursor_mz", ""))))
        if spectrum.get("retention_time"):
            spec_id.append(str(float(spectrum.get("retention_time", ""))))
        spec_id.append(filename)
        spec_id = "_".join(spec_id)
        reference_id = spectrum.get("compound_name", spec_id)
        spectrum = MS2Spectrum(
            spec_id=spec_id,
            precursor_mz=float(spectrum.get("precursor_mz")),
            precursor_rt=float(spectrum.get("retention_time")) if spectrum.get("retention_time") else None,
            list_mz=[],
            list_intensity=[],
            matchms_spectrum=spectrum,
            source=filename,
            instrument="Unknown",
            collision_energy="Unknown",
            compound_name=reference_id,
        )
    else:
        spec_id = []
        spec_id.append(str(float(spectrum["prec_mz"])))
        spec_id.append(str(float(spectrum["rt"])))
        spec_id.append(filename)
        spec_id = "_".join(spec_id)
        precursor_mz = spectrum["prec_mz"]
        precursor_rt = spectrum["rt"]
        reference_id = spectrum["cpd_name"]
        spectrum = spectrum["spectrum"]
        if not skip_filters:
            spectrum = matchms.filtering.default_filters(spectrum)
            spectrum = matchms.filtering.normalize_intensities(spectrum)
            spectrum = matchms.filtering.require_minimum_number_of_peaks(
                spectrum, min_peaks
            )
        spectrum.set("retention_time", precursor_rt)
        spectrum.set("precursor_mz", precursor_mz)

        spectrum = MS2Spectrum(
            spec_id=spec_id,
            precursor_mz=float(precursor_mz),
            precursor_rt=float(precursor_rt),
            list_mz=[],
            list_intensity=[],
            matchms_spectrum=spectrum,
            source=filename,
            instrument="Unknown",
            collision_energy="Unknown",
            compound_name=reference_id,
        )

    return spectrum


def search_for_mzml(sdir):
    """
    Given a directory, encrively search for all files with an .mzml
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
    if isinstance(to_encode, (str, int, float, type(None))):
        return to_encode
    elif isinstance(to_encode, (list, set)):
        return [recursive_encoder(x) for x in to_encode]
    elif isinstance(to_encode, dict):
        return {k: recursive_encoder(v) for k, v in to_encode.items()}
    elif hasattr(to_encode, "json_repr"):
        return to_encode.json_repr
    elif hasattr(to_encode, "serialize"):
        return recursive_encoder(to_encode.serialize())


# the following define mappings between strings and functions used in
# many places in the pcpfm. New modes can be added trivially by adding
# them to these mappings.

correlation_modes = {
    "spearman": scipy.stats.spearmanr,
    "kendall": scipy.stats.kendalltau,
    "pearson": np.corrcoef,
}

descriptive_stat_modes = {
    "median": np.median,
    "max": np.max,
    "min": np.min,
    "mean": np.mean,
}

log_modes = {
    "log2": np.log2, 
    "log10": np.log10
    }

def delete_dir_or_file(to_delete):
    try:
        shutil.rmtree(to_delete)
    except:
        try:
            os.remove(to_delete)
        except:
            "Unable to delete: " + to_delete

file_operations = {
    "link": os.link,
    "copy": shutil.copy,
    "move": shutil.move,
    "delete": delete_dir_or_file
}
