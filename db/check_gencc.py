with open("gencc_submissions.tsv") as f:
    all_ids = {}
    header = []
    for line in f:
        if line.startswith("uuid"):
            header = line.strip().split("\t")
            continue

        columns = line.strip().split("\t")
        line_info = dict(zip(header, columns))

        uuid = line_info["uuid"]
        gene = line_info["gene_curie"].replace(":", "_")
        gene2 = line_info["submitted_as_hgnc_id"].replace(":", "_")
        omim = line_info["disease_original_curie"].replace(":", "_")
        omim2 = line_info["submitted_as_disease_id"].replace(":", "_")

        if uuid.split("-")[1] != gene:
            print(gene, gene2)

        # if new_id in all_ids:
        #     print(all_ids[new_id])
        #     print(uuid)
        # else:
        #     all_ids[new_id] = uuid
