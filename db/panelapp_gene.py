class PanelAppGene:
    def __init__(self, gene):
        self.gene_data = gene["gene_data"]

        self.alias = self.gene_data["alias"]
        self.biotype = self.gene_data["biotype"]
        self.hgnc_id = self.gene_data["hgnc_id"]
        self.gene_name = self.gene_data["gene_name"]
        self.omim_gene = self.gene_data["omim_gene"]
        self.alias_name = self.gene_data["alias_name"]
        self.gene_symbol = self.gene_data["gene_symbol"]
        self.hgnc_symbol = self.gene_data["hgnc_symbol"]
        self.hgnc_release = self.gene_data["hgnc_release"]
        self.ensembl_genes = self.gene_data["ensembl_genes"]

        self.ensembl_genes_grch37 = self.ensembl_genes["GRch37"]
        self.ensembl_genes_grch38 = self.ensembl_genes["GRch38"]
        self.alias = self.gene_data["alias"]
