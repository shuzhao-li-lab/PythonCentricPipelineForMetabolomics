import sys
import pandas as pd
import json

ft = pd.read_csv(sys.argv[1], sep="\t")
empcpd = json.load(open(sys.argv[2]))

annotations = {}
for feature_id in ft["id_number"]:
    annotations[feature_id] = {
        "AnnotFormulas": [],
        "AnnotL4": [],
        "AnnotL2": []
    }

for id, cpd in empcpd.items():
    cpd_formula = ''
    cpd_names = ''
    l2annots = ''
    if "list_matches" in cpd:
        formulas = []
        for x in cpd["list_matches"]:
            formulas.append(x[0].split("_")[0])
        cpd_formula = ",".join(formulas)
    if "mz_only_db_matches" in cpd:
        possible_hits = []
        for x in cpd["mz_only_db_matches"]:
            if x["name"]:
                possible_hits.append(x["name"] + "_" + x["primary_id"])
            else:
                possible_hits.append(x["primary_id"])
        cpd_names = ";".join(possible_hits)
    
    l2s = []
    for x in cpd["MS2_Spectra"]:
        feature_id = x["Matching_Feature"]
        for y in x["Annotations"]:
            msms_score = round(y["msms_score"],3)
            id = y["reference_id"]
            l2s.append(str(msms_score) + "_" + id + "_" + feature_id)
    l2annots = ";".join(l2s)
    for p in cpd["MS1_pseudo_Spectra"]:
        p_id = p["id_number"]
        annotations[p_id]["AnnotFormulas"].append(cpd_formula)
        annotations[p_id]["AnnotL4"].append(cpd_names)
        annotations[p_id]["AnnotL2"].append(l2annots)

AnnotFormulas = []
AnnotL4 = []
AnnotL2 = []
for feature_id in ft["id_number"]:
    AnnotFormulas.append(",".join(annotations[feature_id]["AnnotFormulas"]))
    AnnotL4.append(";".join(annotations[feature_id]["AnnotL4"]))
    AnnotL2.append(";".join(annotations[feature_id]["AnnotL2"]))
ft["AnnotFormulas"] = AnnotFormulas
ft["AnnotL4"] = AnnotL4
ft["AnnotL2"] = AnnotL2

ft.to_csv(sys.argv[1].replace(".tsv", "_annotated.tsv"), sep="\t")
