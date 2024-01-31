from metDataModel.core import Spectrum
import os 

class MS2Spectrum(Spectrum):
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
        source = os.path.basename(source) if source != '' else os.path.basename(source)
        self.precursor_ion = str(precursor_mz) + "_" + str(precursor_rt) + "_" + os.path.basename(source)
        self.retention_time = precursor_rt
        self.precursor_ion_mz = precursor_mz
        self.list_mz = list_mz if list_mz is not None else []
        self.list_intensity = list_intensity if list_intensity is not None else []
        self.instrument = instrument
        self.collision_energy = collision_energy
        self.matchms_spectrum = matchms_spectrum
        if self.matchms_spectrum:
            self.list_mz = [x[0] for x in self.matchms_spectrum.peaks]
            self.list_intensity = [x[1] for x in self.matchms_spectrum.peaks]
        self.annotations = []
        self.compound_name = compound_name
        self.source = source

    def annotate(self, other_MS2, score, matched_peaks, L1_annotation=False):
        self.annotations.append(
            {
                "msms_score": score,
                "matched_peaks": matched_peaks,
                "db_precursor_mz": other_MS2.precursor_ion_mz,
                "reference_id": other_MS2.compound_name,
                "db_spectrum": [[x[0], x[1]] for x in other_MS2.matchms_spectrum.peaks],
                "annot_source": other_MS2.source,
                "annotation_level": "Level_2" if L1_annotation is False else "Level_1"
            }
        )

    def embedding(self):
        embedding = {}
        for k, v in self.__dict__.items():
            if type(v) in [int, float, str, dict, set, list]:
                embedding[k] = v
        return embedding