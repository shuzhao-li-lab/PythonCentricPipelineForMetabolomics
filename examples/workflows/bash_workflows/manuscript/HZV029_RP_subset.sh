working_dir=~/Analyses/
experiment=~/Analyses/HZV029_RP_subset

pcpfm assemble -o $working_dir -j HZV029_RP_subset -s ./subset_sequence_files/RP_neg_subset_sequence.csv --name_field "File Name"
pcpfm asari -i $experiment
pcpfm build_empCpds -i $experiment -tm full -em full
pcpfm generate_output -i $experiment -em full -tm preferred