import matchms

def get_parser(file_extension):
    try:
        return matchms.importing.__getattribute__("load_from_" + file_extension.lower())
    except:
        raise Exception("no matching parser for file type: ", file_extension)

def get_similarity_method(method_name):
    try:
        return matchms.similarity.__getattribute__(method_name)
    except:
        raise Exception("no matching similarity method named: ", method_name)

def process_ms2_spectrum(spectrum, filename="not_specified"):
    min_peaks = 3
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
    

