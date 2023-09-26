# Check PMID from aggrgate data and determine how to preprocess PMID
weirdos = set()
with open("panelapp.uk.aggregate.tsv") as f:
    for line in f:
        if line.startswith("#"):
            continue

        columns = line.strip().split("\t")

        pmids = columns[4].split("||")

        for pmid in pmids:
            if pmid.isnumeric():
                continue

            weirdos.add(pmid)

for weirdo in weirdos:
    print(weirdo)
