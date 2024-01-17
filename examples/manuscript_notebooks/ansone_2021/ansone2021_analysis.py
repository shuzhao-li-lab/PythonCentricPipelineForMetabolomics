import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv 
from matplotlib import colors as mcolors
import json

sys.setrecursionlimit(100000000)

ft = pd.read_csv(sys.argv[1], sep="\t")
ids = ft['id_number']
annots = json.load(open(sys.argv[2]))
feature_to_khipu = {}
for kp_id, kp in annots.items():
    for p in kp["MS1_pseudo_Spectra"]:
        feature_to_khipu[p['id_number']] = kp_id

for column in ft.columns:
    if "Smpl" not in column or "QC"  in column:
        ft.drop(columns=column, inplace=True)
metadata = pd.read_csv("/Users/mitchjo/Downloads/s_MTBLS3852.csv")



colors = ['k', 'red', 'green', 'blue']
color_map = {}
stage_map = {}
color_index = 0
col_colors = []
for column in ft.columns:
    stage = metadata[metadata["Source Name"] == column]["Factor Value[Disease staging]"].values[0]
    if stage not in color_map:
        color_map[stage]=colors[color_index]
        color_index += 1
    col_colors.append(color_map[stage])
    stage_map[column] = stage

def __extract__(row, columns):
    return [row[c] for c in columns]
        

control_cols = [x for x, v in stage_map.items() if v == "control"]
severe_cols = [x for x, v in stage_map.items() if v == "severe COVID-19"]
post_cols = [x for x, v in stage_map.items() if v == "post-acute COVID-19"]


control_vals = ft.apply(__extract__, axis=1, args=(control_cols,))
severe_vals = ft.apply(__extract__, axis=1, args=(severe_cols,))
post_vals = ft.apply(__extract__, axis=1, args=(post_cols,))

from scipy.stats import ttest_ind

p_vals_for_ids = []
for c,s,p, id in zip(control_vals, severe_vals, post_vals, ids):
    t,p = ttest_ind(c, s + p, equal_var=True)
    p_vals_for_ids.append((id, p))

p_vals_for_ids = sorted(p_vals_for_ids, key=lambda x: x[1])
for x in p_vals_for_ids:
    print(x, feature_to_khipu[x[0]])



sns.clustermap(np.log2(ft + 1), col_colors=col_colors)
plt.show()