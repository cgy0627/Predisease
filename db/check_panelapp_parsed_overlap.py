check_overlap = {}
with open("panelapp.uk.genes_filtered.tsv") as f:
    header = []
    for line in f:
        if line.startswith("identifier"):
            header = line.strip().split("\t")
            continue

        columns = line.strip().split("\t")
        line_info = dict(zip(header, columns))

        identifier = line_info["identifier"]
        panel_id = line_info["panel_id"]
        panel_version = line_info["panel_version"]
        tags = line_info["tags"]
        line_number = line_info["line_number"]

        new_id = (
            identifier
            + "-"
            + panel_id
            + "-"
            + panel_version
            + "-"
            + line_number
        )

        if identifier in check_overlap:
            check_overlap[new_id].append(line.strip())
        else:
            check_overlap[new_id] = [line.strip()]

for identifier, line in check_overlap.items():
    if len(line) > 1:
        print(identifier)
        print("\n".join(line))
