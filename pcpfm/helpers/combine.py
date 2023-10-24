import csv
import sys
import os
import multiprocessing as mp
import re

def process(raw, output_dir, existing, delete=True, cleanup=True, add_directory_name=True):
    os.makedirs(output_dir, exist_ok=True)
    try:
        mzML_name = os.path.basename(raw).replace(".raw", ".mzML")
        if add_directory_name:
            mzML_name = os.path.basename(os.path.dirname(raw)) + "___" + mzML_name
        output_file = os.path.join(os.path.abspath(output_dir), mzML_name)
        if mzML_name in existing:
            if cleanup and delete:
                os.remove(raw)
            return True
        else:    
            command = "/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics-1/ThermoRawFileParser/ThermoRawFileParser.exe -f 2 -b " + output_file + " -i " + raw
            try:
                return_val = os.system(command)
                if delete and return_val == 0:
                    os.remove(raw)
            except:
                pass
        return True
    except:
        return False

def find_files(directory):
    queue = [os.path.abspath(directory)]
    raw_files = []
    csv_files = []
    while queue:
        working = queue.pop()
        for file in os.listdir(working):
            abspath = os.path.join(working, file)
            if os.path.isfile(abspath):
                if not file.startswith("._"):
                    if file.endswith(".csv"):
                        csv_files.append(abspath)
                    if file.endswith(".raw"):
                        raw_files.append(abspath)
            elif os.path.isdir(abspath):
                queue.append(abspath)
    return raw_files, csv_files

def combine_csv_files(in_dir, out_dir, csv_out):
    all_entries = [] 
    x = 0
    file_no = 0
    for root, _, files in os.walk(in_dir):
        for file in files:
            filepath = os.path.join(root, file)
            if file.endswith(".csv") and not file.startswith("._"):
                print(file)
                batch_no = re.search("batch_(\d+)", filepath.lower())
                if batch_no is None:
                    batch_no = re.search("batch(\d+)", filepath.lower())
                if batch_no is None:
                    batch_no = "no_batch"
                if batch_no:
                    #batch_no = int(batch_no.group(1))
                    for entry in csv.DictReader(open(filepath).readlines()[1:]):
                        if "File Name" in entry:
                            file_no += 1
                            entry['Filepath'] = os.path.join(os.path.abspath("."), out_dir, entry["File Name"] + ".mzML")
                            if os.path.exists(entry['Filepath']):
                                pass
                            else:
                                path = entry["Path"].split('\\')[-1]
                                #path = "09222023_HMGCS2_Zukai_Cellpellets"
                                entry['Filepath'] = os.path.join(os.path.abspath("."), out_dir, path + "___" + entry["File Name"] + ".mzML")
                            print(entry["Filepath"])
                            if os.path.exists(entry['Filepath']):
                                entry['Method'] = entry['Instrument Method']
                                entry['batch'] = batch_no
                                entry['file_no'] = file_no
                                entry['Sample Type'] = entry["Sample Type"].replace(" ", "_").lower()
                                file_name = entry["File Name"]
                                entry["File Name"] = file_name.lower() + entry["Sample ID"].lower()
                                new_sample_type = []
                                entry['instrument'] = "HFX"
                                if "IDX" in entry["Path"]:
                                    entry['instrument'] = 'IDX'
                                if 'nist' in entry["File Name"] or "nist" in entry["Sample ID"]:
                                    new_sample_type.append("nist")
                                if "qstd" in entry["File Name"]:
                                    new_sample_type.append("qstd")
                                if "pool" in entry["File Name"]:
                                    new_sample_type.append("pooled")
                                if "blank" in entry["File Name"]:
                                    new_sample_type.append("blank")
                                if "sol" in entry["File Name"]:
                                    new_sample_type.append("solvent")
                                if "process" in entry["File Name"]:
                                    new_sample_type.append("process")
                                if "_is" in entry["File Name"] or "std" in entry["File Name"]:
                                    new_sample_type.append("standard")
                                if "dda" in entry["File Name"]:
                                    new_sample_type.append("dda")
                                entry["File Name"] = file_name
                                if new_sample_type:
                                    entry["Sample Type"] = "_".join(new_sample_type)
                                else:
                                    entry["Sample Type"] = "unknown"
                                all_entries.append(entry)
                            else:
                                print(entry["Filepath"])

    all_entries = sorted(all_entries, key=lambda x: (x["batch"], x["file_no"]))
    with open(csv_out, 'w+') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=all_entries[0].keys())
        writer.writeheader()
        for entry in all_entries:
            writer.writerow(entry)
    return all_entries

def main(in_dir, out_dir, csv_out):
    raw_files, _ = find_files(in_dir)
    existing_mzML = set()
    for _, _, files in os.walk(out_dir):
        for file in files:
            if file.endswith(".mzML"):
                existing_mzML.add(file)
    #workers = mp.Pool(1)
    workers = mp.Pool(mp.cpu_count()-1)
    workers.starmap(process, [[r, out_dir, existing_mzML] for r in raw_files])
    _ = combine_csv_files(in_dir, out_dir, csv_out)

if __name__ == '__main__':
    in_dir = sys.argv[1]
    mzML_dir = sys.argv[2]
    csv_out = sys.argv[3]
    main(in_dir, mzML_dir, csv_out)

