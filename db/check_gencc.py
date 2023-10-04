with open("gencc_submissions.tsv") as f:
    all_ids = {}
    prefix = set()
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
        submitter = line_info["submitter_title"]
        submitter_as_submitted = line_info["submitted_as_submitter_name"]

        mondo = line_info["disease_curie"]

        if submitter != submitter_as_submitted:
            print(submitter, submitter_as_submitted)

        # prefix.add(mondo.split(":")[0])
        # if uuid.split("-")[2] != omim:
        #     print(uuid)

        # if new_id in all_ids:
        #     print(all_ids[new_id])
        #     print(uuid)
        # else:
        #     all_ids[new_id] = uuid

for pre in prefix:
    print(pre)
