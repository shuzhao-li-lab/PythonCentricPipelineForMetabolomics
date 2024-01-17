import sys
import pandas as pd
import os
from intervaltree import IntervalTree

standards = pd.read_excel("./CheckMate_annots.xlsx")
pft = None
fft = None
for d, _, fs in os.walk("/Users/mitchjo/Documents/pcpfm_submission/supplemental/Analyses/checkmate_orbi"):
    for f in fs:
        if "preferred" in f:
            pft = pd.read_csv(os.path.join(d, f), sep="\t")
        if "full" in f:
            fft = pd.read_csv(os.path.join(d, f), sep="\t")

standards_tree_mz = IntervalTree()
standards_tree_rt = IntervalTree()
all_names = set()
for name, mz, rtime in zip(standards["Metabolite"], standards["m/z"], standards["Retention time (min)"]):
    mz_error = mz / 1e6 * 10
    standards_tree_mz.addi(mz - mz_error, mz + mz_error, name)
    standards_tree_rt.addi(rtime * 60 - 30, rtime * 60 + 30, name)
    all_names.add(name)

for table in [pft, fft]:
    found = set()
    for mz, rtime in zip(table["mz"], table["rtime"]):
        mz_matches = standards_tree_mz.at(mz)
        rt_matches = standards_tree_rt.at(rtime)
        mz_names = set([x.data for x in mz_matches])
        rt_names = set([x.data for x in rt_matches])
        common = mz_names.intersection(rt_names)
        for n in common:
            found.add(n)
    print(len(found))
    for name in all_names:
        if name not in found:
            print("\tNOT FOUND: ", name)
