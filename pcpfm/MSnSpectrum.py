"""
MS2Spectrum contains the constructors for MS2 spectra. These are used for L1a and L2 annotations. 

The actual MS2 comparisons are performed by matchms in the EmpCpds object.

"""

import os
from matchms.Spectrum import Spectrum
import numpy as np
from metDataModel import core

class MS2Spectrum(core.Spectrum):
    """
    This represents an MS2 spectrum
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
            id (str): a string that uniquely identifies the MS2 spectrum
            precursor_mz (float): the m/z of the precursor ion
            precursor_rt (float): the retention time for the precurso ion
            list_mz (list, optional): list of m/z values for the spectrum. Defaults to None.
            list_intensity (list, optional): list of intensities for the spectrum co-indexed with list_mz. Defaults to None.
            matchms_spectrum (object, optional): a matchms spectrum object. Defaults to None.
            source (str, optional): the file from which the ms2 spectrum was created. Defaults to "".
            instrument (str, optional): the type of instrument from which the spectrum was acquired. Defaults to None.
            collision_energy (float, optional): the collision energey in eV. Defaults to None.
            compound_name (str, optional): the name of the compound. Defaults to None.
            annotations (list, optional): list of annotations for the spectrum. Defaults to None.
        """
        super().__init__(spec_id)
        source = os.path.basename(source) if source != "" else os.path.basename(source)
        self.precursor_ion_id = (str(precursor_mz) + "_" + str(precursor_rt) + "_" + os.path.basename(source))
        self.spec_id = spec_id
        self.rtime = precursor_rt
        self.retention_time = precursor_rt
        self.precursor_ion_mz = precursor_mz
        self.instrument = instrument
        self.collision_energy = collision_energy
        self.matchms_spectrum = matchms_spectrum
        if self.matchms_spectrum:
            self.list_mz = [x[0] for x in self.matchms_spectrum.peaks]
            self.list_intensity = [x[1] for x in self.matchms_spectrum.peaks]
            self.min_mz = min(self.list_mz)
            self.max_mz = max(self.list_mz)
        elif list_mz and list_intensity:
            self.list_mz = list_mz
            self.list_intensity = list_intensity
        else:
            self.list_mz = []
            self.list_intensity = []


        self.annotations = [] if annotations is None else annotations
        self.compound_name = compound_name
        self.source = source

    @property
    def prec_mz(self):
        """
        Simply a shortcut for accessing the precursor_ion_mz field. This is used because the 
        precursor_ion_mz is long and can yield lines that are too long for pylint.

        Returns:
            float: precursor ion mz
        """
        return self.precursor_ion_mz

    @staticmethod
    def from_embedding(embedding):
        """
        This recreates the MS2 object from the serialized version of the object

        Args:
            embedding (dict): the serialized version of the MSnSpectrum.

        Returns:
            MS2Spectrum: the MSnSpectrum that was serialized
        """
        return MS2Spectrum(
            spec_id=None,
            precursor_mz=float(embedding["precursor_ion_mz"]),
            precursor_rt=float(embedding["rtime"]),
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
        """
        Given another MSnSpectrum object, annotate this MSnSpectrum.

        Args:
            other_MS2 (object): the other MS2 spectrum, typically from a database
            score (float): the score from the similarity comparison
            matched_peaks (int): number of peaks 
            annotation_level (str, optional): a string designating the annotation level. Defaults to "Unspecified".
        """
        self.annotations.append(
            {
                "msms_score": score,
                "matched_peaks": matched_peaks,
                "db_precursor_mz": other_ms2.precursor_ion_mz,
                "db_precursor_rt": other_ms2.rtime,
                "reference_id": other_ms2.compound_name,
                "list_mz": [x[0] for x in other_ms2.matchms_spectrum.peaks],
                "list_intensity": [x[1] for x in other_ms2.matchms_spectrum.peaks],
                "primary_db": other_ms2.source,
                "source": self.source,
                "annotation_level": annotation_level,
            }
        )

    def embedding(self):
        """
        Given an MSnSpectrum, generate a serializable version of the object.

        Returns:
            dict: serialized version
        """
        embedding = {}
        for k, v in self.__dict__.items():
            if isinstance(v, (int, float, str, dict, set, list)):
                embedding[k] = v
        return embedding
