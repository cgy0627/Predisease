import sys

input_file = sys.argv[1]

whole_line = ""
with open(input_file) as f:
    for line in f:
        if line.startswith('"uuid'):
            infos = line.strip().split("\t")

            print(infos[0][1:-1], end="")

            for i in range(1, len(infos)):
                if infos[i] == '""':
                    print(f"\t-", end="")
                else:
                    print(f"\t{infos[i][1:-1]}", end="")

            print()

        if not line.startswith('"GENCC'):
            whole_line = whole_line + "\t" + line.strip()
            continue

        if whole_line:
            infos = line.strip().split("\t")

            print(infos[0][1:-1], end="")

            for i in range(1, len(infos)):
                if infos[i] == '""':
                    print(f"\t-", end="")
                else:
                    print(f"\t{infos[i][1:-1]}", end="")

            print()

        whole_line = line.strip()
