import sys
import pymzml
import pandas as pd
from intervaltree import IntervalTree
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity


exp = pymzml.run.Reader(sys.argv[1])
ft = pd.read_csv(sys.argv[2], sep="\t")
ms2_db = open(sys.argv[3])

spectra_registry = {}
for spec in exp:
    if spec.ms_level == 2:
        scan_time = spec.scan_time[0] * 60
        precursor_mz = spec.selected_precursors[0]['mz']
        spectra = spec.peaks('centroided')
        spectra_registry[len(spectra_registry)+1] = (scan_time, precursor_mz, spectra)

rt_tolerance = 20
mz_tolerance = 5

rt_tree = IntervalTree()
mz_tree = IntervalTree()
for id, spec_data in spectra_registry.items():
    scan_time, precursor_mz, spectra = spec_data
    mz_error = precursor_mz * mz_tolerance / 1e6
    rt_tree.addi(max(0, scan_time - rt_tolerance), max(0, scan_time + rt_tolerance), id)
    mz_tree.addi(precursor_mz - mz_error, precursor_mz + mz_error, id)

ms2_matches = []
ms2_to_search = set()
for f_id, f_mz, f_rtime in zip(ft['id_number'], ft['mz'], ft['rtime']):
    mz_matches = mz_tree.at(f_mz)
    if mz_matches:
        rt_matches = rt_tree.at(f_rtime)
        if rt_matches:
            mz_rt_matching_ids = set([x.data for x in mz_matches]).intersection(set([x.data for x in rt_matches]))
            if mz_rt_matching_ids:
                for matched_id in mz_rt_matching_ids:
                    ms2_to_search.add(matched_id)
                ms2_matches.append(mz_rt_matching_ids)
            else:
                ms2_matches.append(None)
        else:
            ms2_matches.append(None)
    else:
        ms2_matches.append(None)

db_entries = []
db_entry = {}
found_peaks = 0
for line in ms2_db:
    if line.startswith("Name: "):
        db_entry["name"] = line.rstrip()
    elif line.startswith("Synon: "):
        db_entry["synon"] = line
    elif line.startswith("PrecursorMZ: "):
        db_entry["precursor_mz"] = float(line.split()[-1])
    elif line.startswith("Num Peaks:"):
        db_entry["num peaks"] = int(line.split()[-1])
        db_entry["spectrum"] = []
    elif line.startswith("Ion_mode: "):
        db_entry["ion_mode"] = line.split()[-1]
    elif "num peaks" in db_entry and found_peaks < db_entry["num peaks"]:
        mz, intensity = line.split()
        found_peaks += 1
        db_entry["spectrum"].append([mz, intensity])
        if found_peaks == db_entry["num peaks"]:
            if "precursor_mz" in db_entry:
                precursor_mz = db_entry["precursor_mz"] 
                if mz_tree.at(precursor_mz):
                    print("match")
                    db_entries.append(db_entry)
                else:
                    print("skip")
            db_entry = {}
            found_peaks = 0

multiplier = 10
exp_spectra = []
for id, spectrum in spectra_registry.items():
    new_spectrum = [0 for _ in range(500 * multiplier)]
    sum_intensity = 0
    for mz, intensity in spectrum[2]:
        if mz < 499:
            sum_intensity += intensity
    for mz, intensity in spectrum[2]:
        if mz < 499:
            new_spectrum[int(round(mz * multiplier,0))] += (intensity / sum_intensity) * 100
    exp_spectra.append(np.array(new_spectrum, dtype=np.float64))
    spectra_registry[id] = (spectrum[0], spectrum[1], new_spectrum)

db_spectra = []
for i, entry in enumerate(db_entries):
    if entry["ion_mode"] == 'P':
        new_spectrum = [0 for _ in range(500 * multiplier)]
        for mz, intensity in entry["spectrum"]:
            if float(mz) < 499:
                sum_intensity += float(intensity)
        for mz, intensity in entry["spectrum"]:
            if float(mz) < 499:
                new_spectrum[int(round(float(mz) * multiplier,0))] += (float(intensity) / sum_intensity) * 100
        db_spectra.append(np.array(new_spectrum, dtype=np.float64))
        entry["spectrum"] = new_spectrum

id_match = {x: 0 for x in range(len(db_spectra))}
for id, exp_spec_data in spectra_registry.items():
    print(id)
    exp_precursor_mz = exp_spec_data[1]
    exp_precursor_mz_error = exp_precursor_mz * 5 / 1e6
    exp_spec = exp_spec_data[2]
    matching_ids = [(x,y) for (x,y) in enumerate(cosine_similarity([exp_spec], db_spectra)[0]) if y > .60]
    if matching_ids:
        for matching_id, score in matching_ids:
            id_match[matching_id] += 1
            found = False
            try:
                db_precursor = db_entries[matching_id]["precursor_mz"]
                if exp_precursor_mz - exp_precursor_mz_error < db_precursor < exp_precursor_mz + exp_precursor_mz_error:
                    print("\t", db_entries[matching_id]["name"], db_entries[matching_id]["precursor_mz"], exp_spec_data[1], score)
                if exp_precursor_mz - exp_precursor_mz_error - 1.003355 < db_precursor < exp_precursor_mz - 1.003355 + exp_precursor_mz_error:
                    print("\t", db_entries[matching_id]["name"], db_entries[matching_id]["precursor_mz"], exp_spec_data[1], score)
            except:
                pass
plt.bar([x for x in sorted(list(id_match.keys()))], [id_match[x] for x in sorted(list(id_match.keys()))])
plt.show()
