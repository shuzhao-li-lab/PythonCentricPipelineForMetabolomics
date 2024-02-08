This directory contains example nextflow workflows for the pcpfm.

Essentiall they are wrappers around bash workflows in `compatible_workflows`` subdirectory. They will copy the results from the generate_output command to the directory where the nextflow script was ran. 

To use these, you need to change pcpfm_workflow param to the absolute path to the script. You will also need to change the jobs field to point to the sequence file you want to use. 

Example jobs are given for the example bash workflows in the `compatible_workflows` directory.

In general, any bash workflow can be converted to a nextflow compatible one by changing the output directory to `.$project_name`. 

