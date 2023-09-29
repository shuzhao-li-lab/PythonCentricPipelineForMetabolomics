pcpfm assemble -s ./new_csv.csv --name_field="File Name" --path_field="InferredPath" -o . -j Test2
pcpfm convert -i ./Test2/
pcpfm asari -i ./Test2/

