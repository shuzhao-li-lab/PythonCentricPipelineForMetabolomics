# these declaration help with processing as the experiment directory can be easily defined
# Working dir needs to be the directory in which to store the analysis directories.
# experiment needs to match the -j for the assemble command

working_dir=~/Analyses
experiment=~/Analyses/checkmate_orbi_subset/

# this simply sources the venv in which the pcpfm is installed.
#source ../pcpfm_venv/bin/activate

# delete any previous analysis
rm -rf $experiment

pcpfm assemble -o $working_dir -j checkmate_orbi_subset -s ./subset_sequence_files/checkmate_orbi_subset_sequence.csv --name_field "File Name"
pcpfm asari -i $experiment
