"""_summary_

"""

import os
from matchms.Spectrum import Spectrum
import numpy as np

class MS2Spectrum:
    """_summary_
    """
    def __init__(
        self,
        spec_id,
        precursor_mz,
        precursor_rt,
        list_mz=None,
        list_intensity=None,
        matchms_spectrum=None,
        source="",
        instrument=None,
        collision_energy=None,
        compound_name=None,
        annotations=None,
    ):
        """_summary_

        Args:
            id (_type_): _description_
            precursor_mz (_type_): _description_
            precursor_rt (_type_): _description_
            list_mz (_type_, optional): _description_. Defaults to None.
            list_intensity (_type_, optional): _description_. Defaults to None.
            matchms_spectrum (_type_, optional): _description_. Defaults to None.
            source (str, optional): _description_. Defaults to "".
            instrument (_type_, optional): _description_. Defaults to None.
            collision_energy (_type_, optional): _description_. Defaults to None.
            compound_name (_type_, optional): _description_. Defaults to None.
            annotations (_type_, optional): _description_. Defaults to None.
        """
        # super().__init__(id)
        source = os.path.basename(source) if source != "" else os.path.basename(source)
        self.precursor_ion = (
            str(precursor_mz) + "_" + str(precursor_rt) + "_" + os.path.basename(source)
        )
        self.spec_id = spec_id
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
        self.annotations = [] if annotations is None else annotations
        self.compound_name = compound_name
        self.source = source

    @staticmethod
    def from_embedding(embedding):
        """_summary_

        Args:
            embedding (_type_): _description_

        Returns:
            _type_: _description_
        """
        return MS2Spectrum(
            spec_id=None,
            precursor_mz=embedding["precursor_ion_mz"],
            precursor_rt=embedding["retention_time"],
            list_mz=embedding["list_mz"],
            list_intensity=embedding["list_intensity"],
            matchms_spectrum=Spectrum(
                mz=np.array(list(embedding["list_mz"])),
                intensities=np.array(embedding["list_intensity"]),
                metadata={},
            ),
            source=embedding["source"],
            instrument=embedding["instrument"],
            collision_energy=embedding["collision_energy"],
            compound_name=None,
            annotations=embedding["annotations"],
        )

    def annotate(self, other_ms2, score, matched_peaks, annotation_level="Unspecified"):
        """_summary_

        Args:
            other_MS2 (_type_): _description_
            score (_type_): _description_
            matched_peaks (_type_): _description_
            annotation_level (str, optional): _description_. Defaults to "Unspecified".
        """
        self.annotations.append(
            {
                "msms_score": score,
                "matched_peaks": matched_peaks,
                "db_precursor_mz": other_ms2.precursor_ion_mz,
                "reference_id": other_ms2.compound_name,
                "list_mz": [x[0] for x in other_ms2.matchms_spectrum.peaks],
                "list_intensity": [x[1] for x in other_ms2.matchms_spectrum.peaks],
                "annot_source": other_ms2.source,
                "annotation_level": annotation_level,
            }
        )

    def embedding(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        embedding = {}
        for k, v in self.__dict__.items():
            if isinstance(v, (int, float, str, dict, set, list)):
                embedding[k] = v
        return embedding
