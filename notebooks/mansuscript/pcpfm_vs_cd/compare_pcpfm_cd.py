import sys
import json
import intervaltree

def read_CD(path):
    all_annots = []
    primary_annot = None
    sub_header = None
    top_header = None
    for line_no, line in enumerate(open(path)):
        line = line.rstrip()
        if line:
            if line_no == 0:
                top_header = line.split(",")
            if line_no > 0:
                if "Checked" in line:
                    if sub_header is None:
                        sub_header = line.split(",")
                else:
                    data = line.split(",")
                    if data[-1] == data[-2] == data[-3]:
                        alt_annot = dict(zip(sub_header, data))
                        primary_annot['alternatives'].append(alt_annot)
                    else:
                        if primary_annot:
                            if primary_annot["MS2"] != "No MS2":
                                all_annots.append(primary_annot)
                        primary_annot = dict(zip(top_header, data))
                        primary_annot['alternatives'] = []
    if primary_annot:
        if primary_annot["MS2"] != "No MS2":
            all_annots.append(primary_annot)
    return all_annots

def read_PCPFM(path):
    annots = json.load(open(path))
    return annots

def map_cd_annots(cd_annots, pcpfm_annots, mz_err=10, rt_err=30):
    feature_mz_tree = intervaltree.IntervalTree()
    feature_rt_tree = intervaltree.IntervalTree()
    feature_to_khipu = {}
    khipu_to_features = {}
    total_features = 0
    for kp_id, kp in pcpcfm_annots.items():
        khipu_to_features[kp_id] = []
        for ms1 in kp["MS1_pseudo_Spectra"]:
            total_features += 1
            id = ms1['id_number']
            feature_to_khipu[id] = kp_id
            khipu_to_features[kp_id].append(id)
            rtime = ms1['rtime']
            mz = ms1['mz']
            mz_min = mz - mz / 1e6 * mz_err
            mz_max = mz + mz / 1e6 * mz_err
            feature_mz_tree.addi(mz_min, mz_max, id)
            feature_rt_tree.addi(rtime - rt_err, rtime + rt_err, id)
    print("total features: ", total_features)
    mapped_annotated = 0
    total_mapped = 0
    features_used = set()
    for cd_annot in cd_annots:
        if cd_annot['Name']:
            try:
                mz = float(cd_annot['m/z'])
                rt = float(cd_annot['RT [min]'])
                rt = rt * 60
                mz_matches = [x.data for x in feature_mz_tree.at(mz)]
                rt_matches = [x.data for x in feature_rt_tree.at(rt)]
                all_match = [x for x in mz_matches if x in rt_matches]
                if all_match:
                    for x in all_match:
                        for z in khipu_to_features[feature_to_khipu[x]]:
                            features_used.add(z)

                    mapped_annotated += 1
                    total_mapped += len(cd_annot["alternatives"]) + 1
            except:
                pass
    print("CD Features: ", len(features_used))
    print("CD Mapped Unique Annotated: ", mapped_annotated)
    print("CD Mapped Total Annotated: ", total_mapped)
    exit()
    mapped_spectra = 0
    total_spectra = 0
    for cd_annot in cd_annots:
        total_spectra += 1
        try:
            mz = float(cd_annot['m/z'])
            rt = float(cd_annot['RT [min]'])
            rt = rt * 60
            mz_matches = [x.data for x in feature_mz_tree.at(mz)]
            rt_matches = [x.data for x in feature_rt_tree.at(rt)]
            all_match = [x for x in mz_matches if x in rt_matches]
            if all_match:
                mapped_spectra += 1
        except:
            pass
    print("CD Total Spectra: ", total_spectra)
    print("CD Mapped Spectra: ", mapped_spectra)

def count_pcpfm_annots(pcpfm_annots):
    total_mapped_ms2 = 0
    top_mapped_annots = 0
    all_mapped_annots = 0
    features_annotated = set()
    for kp_id, kp in pcpfm_annots.items():
        if "MS2_Spectra" in kp:
            for ms2 in kp["MS2_Spectra"]:
                total_mapped_ms2 += 1
                if ms2["annotations"]:
                    top_mapped_annots += 1
                    all_mapped_annots += min(len(ms2["annotations"]), 10)
                    for peak in kp["MS1_pseudo_Spectra"]:
                        features_annotated.add(peak['id_number'])
    print("PCPFM Asari Annotated: ", len(features_annotated))

    print("PCPFM Total Mapped: ", total_mapped_ms2)
    print("PCPFM Unique Annot Mapped: ", top_mapped_annots)
    print("PCPFM Total Annot Mapped: ", all_mapped_annots)

def count_cd_annots(cd_annots):
    annots_unique = 0
    annots_total = 0
    for cd_annot in cd_annots:
        if cd_annot['Name']:
            annots_unique += 1
            annots_total += 1 + len(cd_annot["alternatives"])
    print("CD Unique Annots: ", annots_unique)
    print("CD Total Annots: ", annots_total)

def compare_annots(cd_annots, pcpfm_annots):
    best_annots_pcpfm = []
    all_annots_pcpfm = []
    for kp_id, kp in pcpfm_annots.items():
        annot_ranked = []
        if "MS2_Spectra" in kp:
            for ms2 in kp["MS2_Spectra"]:
                if ms2["annotations"]:
                    for annot in ms2["annotations"]:
                        annot_ranked.append((annot['reference_id'], annot['msms_score'], ms2['retention_time']))
                    annot_ranked = sorted(annot_ranked, key=lambda x: -x[1])
        if annot_ranked:
            best_annots_pcpfm.append(annot_ranked[0])
            for i in range(min(10, len(annot_ranked))):
                all_annots_pcpfm.append(annot_ranked[i])
                
    
    best_annots_cd = []
    for cd_annot in cd_annots:
        if cd_annot["Name"]:
            try:
                best_annots_cd.append((cd_annot["Name"], float(cd_annot["mzVault Best Match"])/100, 60*float(cd_annot['RT [min]'])))
            except:
                pass
    best_annots_cd = set(best_annots_cd)
    best_annots_pcpfm = set(best_annots_pcpfm)
    all_annots_pcpfm = set(all_annots_pcpfm)

    from sklearn.metrics import r2_score
    cd_common = 0
    Xs, Ys = [], []
    for cd_name, cd_score, cd_time in best_annots_cd:
        best_match = (None, None, 99999999)
        for pcpfm_name, pcpfm_score, pcpfm_time in best_annots_pcpfm:
            if pcpfm_name == cd_name:
                if abs(cd_time - pcpfm_time) < best_match[2]:
                    best_match = (pcpfm_name, pcpfm_score, abs(cd_time - pcpfm_time))
        if best_match[0]: 
            Xs.append(cd_score)
            Ys.append(best_match[1])
            cd_common += 1
    import scipy
    import matplotlib.pyplot as plt
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Xs, Ys)
    print(r_value)
    plt.scatter(Xs, Ys)
    plt.xlabel("CD Score")
    plt.ylabel("PCPFM Score")
    plt.show()

    print("CD Common (top): ", cd_common, len(best_annots_cd))

    cd_common = 0
    Xs, Ys = [], []
    for cd_name, cd_score, cd_time in best_annots_cd:
        best_match = (None, None, 99999999)
        for pcpfm_name, pcpfm_score, pcpfm_time in all_annots_pcpfm:
            if pcpfm_name == cd_name:
                if abs(cd_time - pcpfm_time) < best_match[2]:
                    best_match = (pcpfm_name, pcpfm_score, abs(cd_time - pcpfm_time))
        if best_match[0]: 
            Xs.append(cd_score)
            Ys.append(best_match[1])
            cd_common += 1
    import matplotlib.pyplot as plt
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Xs, Ys)
    print(r_value)
    plt.scatter(Xs, Ys)
    plt.show()
    print("CD Common (top): ", cd_common, len(best_annots_cd))





    

if __name__ == '__main__':
    cd_annots = read_CD(sys.argv[1])
    pcpcfm_annots = read_PCPFM(sys.argv[2])
    count_cd_annots(cd_annots)
    count_pcpfm_annots(pcpcfm_annots)
    compare_annots(cd_annots, pcpcfm_annots)
    map_cd_annots(cd_annots, pcpcfm_annots)