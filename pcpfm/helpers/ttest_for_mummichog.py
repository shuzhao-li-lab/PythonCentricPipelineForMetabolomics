import numpy as np
import stats
import pandas as pd
import sys

def ttest(row, group_a, group_b, name_field="File Name"):
    values_A = np.log2([float(x+1) for x in row[[x[name_field] for x in group_a]] if not np.isnan(x)])
    values_B = np.log2([float(x+1) for x in row[[x[name_field] for x in group_b]] if not np.isnan(x)])
    t, p = stats.ttest_ind(values_A, values_B, equal_var=False)
    if np.isnan(p):
        t, p = 0, 1
    return t, p

def extract(row, columns):
    return {column: row[column] for column in columns}

def to_groups(entries, metadata_field):
    groups = {}
    for entry in entries:
        value = entry[metadata_field]
        if value not in groups:
            groups[value] = []
        groups[value].append(entry)
    return groups


ft = pd.read_csv(sys.argv[1], sep="\t")
metadata = pd.read_csv(sys.argv[2], sep="\t")
metadata_field = sys.argv[3]

metadata_entries = metadata.apply(extract, axis=1, args=(metadata.columns,))
metadata_groups = to_groups(metadata_entries)

used = set()
for g1 in metadata_groups.keys():
    for g2 in metadata_groups.keys():
        if g1 != g2 and (g1, g2) not in used and (g2, g1) not in used:
            for_mummichog = pd.DataFrame()
            used.add((g1, g2))
            used.add((g2, g1))
            for_mummichog["mz"] = ft["mz"]
            for_mummichog["rtime"] = ft["rtime"]
            ttest_results = for_mummichog.apply(ttest, axis=1, args=(g1, g2))
            for_mummichog["t-score"] = [x[0] for x in ttest_results]
            for_mummichog["p-value"] = [x[1] for x in ttest_results]
            for_mummichog.to_csv(g1 + "_" + g2 + ".tsv", sep="\t")