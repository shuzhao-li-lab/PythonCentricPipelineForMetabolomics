from matchms.importing import load_from_mgf, load_from_mzml
from matchms.filtering import default_filters, normalize_intensities, add_precursor_mz
from matchms.filtering import require_minimum_number_of_peaks
from matchms.similarity import CosineGreedy
from matchms.importing import load_from_msp
import pandas as pd
import sys
from intervaltree import IntervalTree
import json

RT_TOL = 5
MZ_TOL = 5
match_MZ = True
match_RT = False
skip_IN_SOURCE=False
SIMILARITY_THRESHOLD=0
MATCH_PEAK_THRESHOLD=3

precursor_mz_tree = IntervalTree()
precursor_rt_tree = IntervalTree()
spectral_registry = {}

file = list(load_from_mzml(sys.argv[2], metadata_harmonization=False))
spectrums = []
for spectrum in file:
    spectrum.set('retention_time', spectrum.metadata['scan_start_time'][0] * 60)
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = add_precursor_mz(spectrum)
    spectrum = require_minimum_number_of_peaks(spectrum, n_required=3)
    if spectrum:
      try:
        if spectrum.get('precursor_mz'):
          precursor_mz = float(spectrum.get('precursor_mz'))
        else:
          precursor_mz = None
        if spectrum.get('retention_time'):
          precursor_rt = float(spectrum.get('retention_time'))
        else:
          precursor_rt = None
        precursor_mz_error = precursor_mz / 1e6 * MZ_TOL
        entry_id = len(spectral_registry)
        spectral_registry[entry_id] = spectrum
        precursor_mz_tree.addi(precursor_mz - precursor_mz_error, precursor_mz + precursor_mz_error, entry_id)
        precursor_rt_tree.addi(precursor_rt - RT_TOL, precursor_rt + RT_TOL, entry_id)
      except:
        pass

#for k, v in spectral_registry.items():
#  print(k, " : ", v)

cosinegreedy = CosineGreedy()
MSMS_annots = []
file_msp = sys.argv[1]

for db_spectrum in load_from_msp(file_msp, metadata_harmonization=False):
  if db_spectrum is not None:
      precursor_mz = None
      if db_spectrum.get('precursor_mz') is not None:
        try:
          precursor_mz = float(db_spectrum.get('precursor_mz'))
        except:
          precursor_mz = None
      #if db_spectrum.get('retention_time') is not None and db_spectrum.get('retention_time'):
      #  print(db_spectrum.metadata)
      #  precursor_rt = float(spectrum.get('retention_time')[0])
      #else:
      #  precursor_rt = None
      mz_matches = set()
      rt_matches = set()
      if match_MZ and precursor_mz:
        for match in precursor_mz_tree.at(precursor_mz):
          mz_matches.add(match.data)
      if match_RT and precursor_rt:
        for match in precursor_rt_tree.at(precursor_rt_tree):
          rt_matches.add(match.data)
      all_matches = set()
      if match_MZ and match_RT:
        all_matches = mz_matches.intersection(rt_matches)
      elif match_MZ:
        all_matches = mz_matches
      elif match_RT:
        all_matches = rt_matches
      else:
        all_matches = set(list(spectral_registry.keys()))

      for experimental_spectrum in [spectral_registry[x] for x in all_matches]:
        msms_score, n_matches = cosinegreedy.pair(db_spectrum, experimental_spectrum).tolist()
        if msms_score > SIMILARITY_THRESHOLD and n_matches > MATCH_PEAK_THRESHOLD:
          feature_mz = experimental_spectrum.get("precursor_mz") if experimental_spectrum.get("precursor_mz") else None
          feature_rt = experimental_spectrum.get("retention_time") if experimental_spectrum.get("retention_time") else None
          match_mz = db_spectrum.get("precursor_mz") if db_spectrum.get("precursor_mz") else None
          match_rt = db_spectrum.get("retention_time") if db_spectrum.get("retention_time") else None

          try:
            match_mz = float(match_mz)
          except:
            match_mz = None

          MSMS_annots.append({
            'msms_score':msms_score,
            'matched_peaks':n_matches,
            'feature_id': None,
            'feature_precursor_mz': feature_mz,
            'feature_precursor_rt': feature_rt,
            'match_precursor_mz': match_mz,
            'match_precursor_rt': match_rt,
            'reference_id': db_spectrum.get("compound_name")})
          print(json.dumps(MSMS_annots[-1], indent=4))
df = pd.DataFrame(MSMS_annots)
df.to_csv('result.csv', sep=',')

