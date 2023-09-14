no_omim = set()
with open("omim.tsv") as omim:
    for line in omim:
        if line.startswith("#"):
            continue

        info = line.strip().split("\t")

        mim = info[0]
        prefix = info[1]
        gene = info[3]
        omim = info[16]
        pmid = info[23]

        if omim == "-":
            no_omim.add(gene)

with open("gencc.no_omim.tsv") as gencc:
    header = []
    for line in gencc:
        if line.startswith("uuid"):
            header = line.strip().split("\t")
            continue

        columns = line.strip().split("\t")
        line_info = dict(zip(header, columns))

        uuid = line_info["uuid"]
        gene = line_info["gene_symbol"]
        mondo_id = line_info["disease_curie"]  # MONDO id with :
        classification = line_info["classification_title"]
        pmid = line_info["submitted_as_pmids"]
        report = line_info["submitted_as_public_report_url"]
        submitter = line_info["submitter_title"]

        hgnc_id = uuid.split("-")[1].strip()  # underscore
        disease_id = uuid.split("-")[2].strip()  # underscore

        identifier = hgnc_id + "-" + disease_id

        # if omim != '-' and gene in no_omim:
        if classification == "Definitive" and gene in no_omim:
            print(
                gene,
                classification,
                pmid,
                disease_id,
                report,
                submitter,
                sep="\t",
            )
