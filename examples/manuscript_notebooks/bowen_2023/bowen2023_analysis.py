import pandas as pd
import sys
import json
from intervaltree import IntervalTree

ft = pd.read_csv(sys.argv[1], sep="\t")
annots = json.load(open(sys.argv[2]))

feature_to_khipu = {}
for kp_id, kp in annots.items():
    for p in kp["MS1_pseudo_Spectra"]:
        feature_to_khipu[p['id_number']] = kp_id
print(feature_to_khipu)

expected_features_cell = [
    (399.1823, 33.4,  'M0'),
    (399.2184, 224.8, 'M0'),
    (371.1874, 257.2, 'M1'),
    (415.2134, 264.8, 'M2'),
    (415.2128, 224.5, 'M2'),
    (343.1562, 249.4, 'M3'),
    (387.1823, 35.4,  'M4'),
    (387.1824, 242.5, 'M4'),
    (413.1978, 34.6,  'M12'),
    (385.1667, 38.9,  'M14'),
    (385.1666, 65.8,  'M14'),
    (159.1490, 273.2, 'M20')
]

expected_features_media = [
    (399.2190, 134.9, 'M0'),
    (399.2191, 224.7, 'M0'),
    (371.1879, 257.1, 'M1'),
    (415.2130, 197.7, 'M2'),
    (415.2141, 264.8, 'M2'),
    (413.1983, 34.5,  'M12'),
    (385.1671, 39.0,  'M14'),
    (344.1407, 35.1,  'M20')
]

mz_tree = IntervalTree()
rt_tree = IntervalTree()

ppm_mz_tol = 5
rt_tol = 5

for mz, rt, id in zip(ft['mz'], ft['rtime'], ft['id_number']):
    mz_err = mz/1e6 * ppm_mz_tol
    mz_tree.addi(mz-mz_err, mz+mz_err, id)
    rt_tree.addi(rt-rt_tol, rt+rt_tol, id)

print("Cell")
for expected in expected_features_cell:
    exp_mz, exp_rt, id = expected
    matches_mz = set([x.data for x in mz_tree.at(exp_mz)])
    matches_rt = set([x.data for x in rt_tree.at(exp_rt)])
    true_matches = matches_mz.intersection(matches_rt)
    print(id, true_matches)

print("\n----\n")
print("Media")
for expected in expected_features_media:
    exp_mz, exp_rt, id = expected
    matches_mz = set([x.data for x in mz_tree.at(exp_mz)])
    matches_rt = set([x.data for x in rt_tree.at(exp_rt)])
    true_matches = matches_mz.intersection(matches_rt)
    print(id, true_matches)