no_omim = set()
with open('omim.tsv') as omim:
    for line in omim:
        if line.startswith("#"):
            continue
        
        info = line.strip().split('\t')

        mim = info[0]
        prefix = info[1]
        gene = info[3]
        omim = info[16]
        pmid = info[23]

        if omim == '-':
            no_omim.add(gene)

with open("panelapp.uk.aggregate.highest_panels.tsv") as f:
    for line in f:
        if line.startswith("#"):
            continue
        
        info = line.strip().split('\t')

        gene = info[1]
        level = info[2]
        pmid = info[4]
        omim = info[5]
        panel_id = info[7]
        panel_version = info[9]

        # if omim != '-' and gene in no_omim:
        if level == '3' and gene in no_omim:
            print(gene, level, pmid, omim, panel_id, panel_version, sep='\t')