#!/bin/bash

# This is the minimal processing example where asari is ran and empcpds are produced but that is it. 

process() {
    pcpfm assemble -o . -j $2 -s $1 --path_field Filepath --name_field "File Name"
    pcpfm asari -i ./$2 
    pcpfm build_empCpds -i ./$2 -tm full -em full 
    pcpfm generate_output -i ./$2 -em full -tm preferred 
    pcpfm report -i ./$2 --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
}

# $1 = sequence file
# $2 = project name
process $1 $2