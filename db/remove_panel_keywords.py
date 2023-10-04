count = 0
output = open("panel_keyword_removed.tsv", "w")
with open("panelapp.uk.aggregate.highest_panels.tsv") as f:
    for line in f:
        if line.startswith("#"):
            print(line.strip(), file=output)
            continue

        info = line.strip().split("\t")

        disease_name = info[8]

        if (
            disease_name == "COVID-19 research"
            or "tumor" in disease_name
            or "tumour" in disease_name
        ):
            count += 1
            continue

        print(line.strip(), file=output)

print(count)
