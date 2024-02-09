working_dir=~/Analyses/
experiment=~/Analyses/HZV029_Lipidomics/

source ../pcpfm_venv/bin/activate

#rm -rf $experiment

pcpfm assemble -o $working_dir -j HZV029_Lipidomics -s ~/Datasets/ForPCPFM/MerckLipidomics_Batch1/HILICneg_RPpos/HILICneg_RPpos.csv --filter /Users/mitchjo/Projects/PythonCentricPipelineForMetabolomics/pcpfm/filters/hilicneg.json --name_field "File Name"
pcpfm convert -i $experiment
pcpfm asari -i $experiment

pcpfm blank_masking --table_moniker preferred --new_moniker blank_masked1 --blank_value Process_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm blank_masking --table_moniker blank_masked1 --new_moniker blank_masked2 --blank_value Solvent_Blank --sample_value Unknown --query_field "Sample Type" --blank_intensity_ratio 3 -i $experiment 
pcpfm drop_samples --table_moniker blank_masked2 --new_moniker masked_preferred_unknowns --drop_value Unknown --drop_field "Sample Type" --drop_others true -i $experiment 
pcpfm drop_outliers --table_moniker masked_preferred_unknowns --new_moniker qaqc_filtered_masked_pref_unknowns -i $experiment
