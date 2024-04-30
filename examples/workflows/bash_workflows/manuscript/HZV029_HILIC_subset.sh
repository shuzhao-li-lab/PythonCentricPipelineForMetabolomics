working_dir=~/Analyses/
experiment=~/Analyses/HZV029_HILIC_subset

pcpfm assemble -o $working_dir -j HZV029_HILIC_subset -s ./subset_sequence_files/HILIC_pos_subset_sequence.csv --name_field "File Name"
pcpfm asari -i $experiment
pcpfm build_empCpds -i $experiment -tm full -em full
pcpfm generate_output -i $experiment -em full -tm preferred