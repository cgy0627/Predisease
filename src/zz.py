with open("pheno") as f:
    for line in f:
        info = line.strip().split("\t")
        phenos = set(info[2].split("||"))

        # if len(phenos) > 1:
        #     print(line.strip())

        pheno = phenos.pop()
        if info[0] != pheno and pheno != "-":
            print(info[0], pheno)
            print(line.strip())
