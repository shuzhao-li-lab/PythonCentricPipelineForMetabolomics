# these declaration help with processing as the experiment directory can be easily defined
# Working dir needs to be the directory in which to store the analysis directories.
# experiment needs to match the -j for the assemble command

working_dir=~/Analyses
experiment=~/Analyses/HZV029_QC/

# this simply sources the venv in which the pcpfm is installed.
#source ../pcpfm_venv/bin/activate

# delete any previous analysis
rm -rf $experiment

pcpfm assemble -o $working_dir -j HZV029_QC -s ./sequence_files/HZV029_QC_sequence.csv --name_field "File Name"
pcpfm convert -i $experiment
pcpfm asari -i $experiment

pcpfm blank_masking --table_moniker preferred --new_moniker blank_masked1 --blank_value Process_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm blank_masking --table_moniker blank_masked1 --new_moniker blank_masked2 --blank_value Solvent_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker blank_masked2 --new_moniker masked_preferred_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
