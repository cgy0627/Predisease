all_ids = {}
with open("panelapp.uk.aggregate.tsv") as f:
    header = []
    for line in f:
        if line.startswith("#"):
            header = line.strip().replace("#", "").split("\t")
            continue

        columns = line.strip().split("\t")
        line_info = dict(zip(header, columns))

        hgnc_id = line_info["hgnc_id"]
        panel_id = line_info["panel_id"]
        panel_version = line_info["panel_version"]

        key = hgnc_id + "-" + panel_id

        if key in all_ids:
            print(key, panel_version, all_ids[key], sep="\t")
        else:
            all_ids[key] = panel_version
