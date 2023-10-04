import sys, os

input_file_path = sys.argv[1]
start_of_line = sys.argv[2]

output_file_base = os.path.basename(input_file_path).split(".")[0]
output_file_path = output_file_base + ".col2idx.txt"

output = open(output_file_path, "w")

with open(input_file_path) as f:
    for line in f:
        if line.startswith(start_of_line):
            header = line.strip().split("\t")

            for i, head in enumerate(header):
                print(i, head, sep="\t", file=output)

            break

output.close()
