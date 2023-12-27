from metDataModel.core import MSnSpectrum
import os 

class MS2Spectrum(MSnSpectrum):
    def __init__(self,
                 id,
                 precursor_mz, 
                 precursor_rt,
                 list_mz=None,
                 list_intensity=None,
                 matchms_spectrum=None,
                 source='',
                 instrument=None,
                 collision_energy=None,
                 compound_name=None):
        super().__init__(id)
        source = os.path.basename(source) if source != '' else source
        self.precursor_ion = str(precursor_mz) + "_" + str(precursor_rt) + "_" + source
        self.retention_time = precursor_rt
        self.precursor_ion_mz = precursor_mz
        self.list_mz = list_mz if list_mz is not None else []
        self.list_intensity = list_intensity if list_intensity is not None else []
        self.instrument = instrument
        self.collision_energy = collision_energy
        self.matchms_spectrum = matchms_spectrum
        self.annotations = []
        self.compound_name = compound_name
        self.source = source

    def annotate(self, other_MS2, score, matched_peaks):
        self.annotations.append(
            {
                "msms_score": score,
                "matched_peaks": matched_peaks,
                "db_precursor_mz": other_MS2.precursor_ion_mz,
                "reference_id": other_MS2.compound_name,
                "db_spectrum": [[x[0], x[1]] for x in other_MS2.matchms_spectrum.peaks],
                "annot_source": other_MS2.source
            }
        )

    def embedding(self):
        embedding = {}
        for k, v in self.__dict__.items():
            if type(v) in [int, float, str, dict, set]:
                embedding[k] = v
        return embedding