#!/bin/bash

# This is the minimal processing example where asari is ran and empcpds are produced but that is it. 

process() {
    experiment=$(pwd)/$2
    pcpfm assemble -o $(pwd) -j $2 -s $1  --path_field Filepath --name_field "File Name"
    pcpfm asari -i $experiment 
    pcpfm build_empCpds -i $experiment -tm full -em full 
    pcpfm generate_output -i $experiment -em full -tm preferred 
    pcpfm report -i $experiment --color_by='["Sample Type"]' --marker_by='["Sample Type"]'
}

# $1 = sequence file
# $2 = project name
process $1 $2