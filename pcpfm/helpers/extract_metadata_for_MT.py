import pandas as pd
import sys

try:
    seq = pd.read_csv(sys.argv[1])
    seq["Sample ID"]
except:
    with open("tmp", 'w+') as temp:
        with open(sys.argv[1]) as sequence:
            for line in sequence.readlines()[1:]:
                temp.write(line)
    seq = pd.read_csv("tmp")

group = []
for x in seq["Sample ID"]:
    group_found = False
    for y in ["CTRL_CPA_MET", "CTRL_CPA", "CTRL", "KD_CPA_MET", "KD_CPA", "KD_MET"]:
        if y in x:
            group_found = True
            group.append(y)
            break
    if not group_found:
        group.append(None)
seq["Group"] = group
seq.to_csv(sys.argv[1].replace(".csv", "_grouped.csv"))