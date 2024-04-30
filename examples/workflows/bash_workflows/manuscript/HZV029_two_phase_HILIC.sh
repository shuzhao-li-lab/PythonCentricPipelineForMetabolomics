working_dir=~/Analyses/
experiment=~/Analyses/HZV029_TwoPhase_HILICneg/

source ../pcpfm_venv/bin/activate

#rm -rf $experiment

pcpfm assemble -o $working_dir -j HZV029_TwoPhase_HILICneg -s ./sequence_files/HZV029_TwoPhase.csv --filter ../../../filters/hilicneg.json --name_field "File Name"
pcpfm convert -i $experiment
pcpfm asari -i $experiment

pcpfm blank_masking --table_moniker preferred --new_moniker blank_masked1 --blank_value Process_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm blank_masking --table_moniker blank_masked1 --new_moniker blank_masked2 --blank_value Solvent_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker blank_masked2 --new_moniker masked_preferred_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
