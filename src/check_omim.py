phenotypes = {}
with open("omim.tsv") as omim:
    header = []
    for line in omim:
        if line.startswith("#"):
            header = line[1:].strip().split("\t")
            continue

        columns = line.strip().split("\t")

        line_info = dict(zip(header, columns))

        mim_number = line_info["mimNumber"]
        prefix = line_info["prefix"]
        gene_symbol = line_info["approvedGeneSymbols"]

        pheno_mim = line_info["genePhenoMap:phenotypeMimNumber"]
        inheritance = line_info["genePhenoMap:phenoInheritance"]
        pubmed_ids = line_info["pubmedId"]

        if prefix == "%":
            phenotypes[mim_number] = 1

        if prefix == "%" and pheno_mim == "-" and pubmed_ids == "-":
            print(mim_number, gene_symbol, pheno_mim)

# with open("omim.tsv") as omim:
#     for line in omim:
#         if line.startswith("#"):
#             header = line[1:].strip().split("\t")
#             continue

#         columns = line.strip().split("\t")

#         line_info = dict(zip(header, columns))

#         mim_number = line_info["mimNumber"]
#         prefix = line_info["prefix"]
#         gene_symbol = line_info["approvedGeneSymbols"]

#         pheno_mim = line_info["genePhenoMap:phenotypeMimNumber"]
#         pubmed_ids = line_info["pubmedId"]

#         if prefix == "" or pheno_mim == "-":
#             continue

#         for mim in pheno_mim.split("||"):
#             if mim in phenotypes:
#                 print(line.strip())
