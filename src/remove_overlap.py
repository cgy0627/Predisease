import requests, yaml, json
from panelapp_gene import PanelAppGene
import re

config_file = "config_temp.yaml"
with open(config_file, "r") as file:
    config = yaml.safe_load(file)


# def parse_mondo_json() -> dict:
#     mondo_filepath = config["database"]["MONDO_JSON_PATH"]
#     mondo_open = open(mondo_filepath, "r")

#     mondo = json.load(mondo_open)
#     # nodes = mondo[0]["graphs"]["nodes"]

#     print(len(mondo["graphs"][0]["nodes"]))


def check_if_omim_id(omim_id: str) -> bool:
    mim = omim_id.split(":")[1]

    return mim.isnumeric()


def parse_mondo() -> dict:
    mondo_filepath = config["database"]["MONDO_PATH"]

    mondo_match = {}
    with open(mondo_filepath, "r") as mondo:
        mondo_id = ""
        omim_ids = set()
        for line in mondo:
            if line.startswith("id: "):
                if omim_ids:
                    mondo_match[mondo_id] = omim_ids

                mondo_id = line.strip().split()[1]
                omim_ids = set()

                if not mondo_id.startswith("MONDO"):
                    mondo_id = ""
                    continue

            if mondo_id == "":
                continue

            is_equivalent = False

            # synonym line
            if line.startswith("synonym: "):
                synonym_line_info = line.strip().split()

                for info in synonym_line_info:
                    info = info.replace("[", "").replace("]", "")

                    if info == "EXACT":
                        is_equivalent = True

                    if (
                        is_equivalent
                        and info.startswith("OMIM:")
                        and check_if_omim_id(info)
                    ):
                        omim_ids.add(info)

            # xref line
            if line.startswith("xref: "):
                xref_line_info = line.strip().split()

                xref_id = xref_line_info[1]

                if len(xref_line_info) < 3:
                    continue

                additional_omim_id = ""
                for additional_info in xref_line_info[2:]:
                    additional_info += (
                        "="  # for cases where = doesn't already exist
                    )

                    source = additional_info.split("=")[1]
                    source = (
                        source.replace('"', "")
                        .replace(",", "")
                        .replace("}", "")
                    )

                    if source == "MONDO:equivalentTo":
                        is_equivalent = True

                    if source.startswith("OMIM:") and check_if_omim_id(source):
                        additional_omim_id = source

                if is_equivalent:
                    if xref_id.startswith("OMIM:"):
                        omim_ids.add(xref_id)
                    if additional_omim_id:
                        omim_ids.add(additional_omim_id)

    # output = open("mondo_match.txt", "w")
    # for id, omims in mondo_match.items():
    #     print(id, omims, file=output)
    # output.close()

    return mondo_match


def parse_omim(omim_db: str) -> dict:
    """Parse OMIM database and save necessary information as a dictionary.

    Note:
        Requires OMIM database given by Sean (omim.tsv).

    Args:
        omim_db (str): Full path for OMIM database.

    Returns:
        dict: Phenotype information for each entry in OMIM database.

    Examples:
        >>> omim_data
        {
            "142622":{
                "prefix":"*",
                "gene_symbol":"HPCA",
                "hpo_strs":"--",
                "phenomap_phenotype":"224500"
            },
            "300905":{
                "prefix":"#",
                "gene_symbol":"-",
                "hpo_strs":"One family has been reported||Decreased auditory brainstem responses",
                "phenomap_phenotype":"300905"
            }
        }
    """
    omim_data = {}
    with open(omim_db) as db:
        header = []
        for line in db:
            if line.startswith("#"):
                header = line.replace("#", "").strip().split("\t")
                continue

            line_info = dict(zip(header, line.strip().split("\t")))

            mim_number = line_info["mimNumber"]

            prefix = line_info["prefix"]
            gene_symbol = line_info["approvedGeneSymbols"]
            hgnc_id = line_info["hgncId"]
            phenomap_phenotype = line_info["genePhenoMap:phenotypeMimNumber"]
            pubmed_id = line_info["pubmedId"]

            if prefix == "-" or prefix == "^":
                continue

            omim_data[mim_number] = {}

            omim_data[mim_number]["prefix"] = prefix
            omim_data[mim_number]["gene_symbol"] = gene_symbol
            omim_data[mim_number]["hgnc_id"] = hgnc_id
            omim_data[mim_number]["phenomap_phenotype"] = phenomap_phenotype
            omim_data[mim_number]["pubmed_id"] = pubmed_id

    return omim_data


def get_omim_ids():
    omim_file_path = "/home/ec2-user/project/Predisease/db/omim.tsv"

    omim_data = parse_omim(omim_file_path)
    omim_ids = set()

    for omim_info in omim_data.values():
        hgnc_id = omim_info["hgnc_id"]
        pheno_mims = omim_info["phenomap_phenotype"]
        pubmed_id = omim_info["pubmed_id"]

        if pheno_mims == "-":
            continue

        for pheno_mim in pheno_mims.split("||"):
            omim_id = "HGNC_" + hgnc_id + "-" + "OMIM_" + pheno_mim

            omim_ids.add(omim_id)

        for pmid_number in pubmed_id.split("||"):
            omim_id = "HGNC_" + hgnc_id + "-" + "PMID_" + pmid_number

            omim_ids.add(omim_id)

    return omim_ids


def remove_omim_from_gencc():
    omim_ids = get_omim_ids()
    mondo_to_omim = parse_mondo()

    gencc_file_path = config["database"]["GENCC_PATH"]

    gencc_no_omim_filepath = config["database"]["GENCC_NO_OMIM_PATH"]
    gencc_no_omim_output = open(gencc_no_omim_filepath, "w")

    gencc_omim_overlap_filepath = "gencc.omim_overlap.tsv"
    gencc_omim_overlap_output = open(gencc_omim_overlap_filepath, "w")

    with open(gencc_file_path) as gencc:
        header = []
        for line in gencc:
            if line.startswith("uuid"):
                header = line.strip().split("\t")
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            uuid = line_info["uuid"]
            mondo_id = line_info["disease_curie"]  # MONDO id with :

            hgnc_id = uuid.split("-")[1].strip()  # underscore
            disease_id = uuid.split("-")[2].strip()  # underscore

            identifier = hgnc_id + "-" + disease_id

            is_novel = True
            if disease_id.startswith("OMIM"):
                if identifier in omim_ids:
                    is_novel = False

            else:
                if mondo_id not in mondo_to_omim:
                    continue

                for temp_omim_id in mondo_to_omim[mondo_id]:
                    identifier_variation = (
                        hgnc_id + "-" + temp_omim_id.replace(":", "_")
                    )

                    if identifier_variation in omim_ids:
                        is_novel = False
                        break

            if is_novel:
                print(line.strip(), file=gencc_no_omim_output)
            else:
                print(line.strip(), file=gencc_omim_overlap_output)

    gencc_no_omim_output.close()
    gencc_omim_overlap_output.close()

    return


def parse_panelapp_aus_json():
    panelapp_output_filepath = config["database"]["PANELAPP_AUS_PATH"]
    panelapp_output = open(panelapp_output_filepath, "w")

    print(
        "#alias",
        "biotype",
        "hgnc_id",
        "gene_name",
        "omim_gene",
        "alias_name",
        "gene_symbol",
        "hgnc_symbol",
        "hgnc_release",
        "ensembl_genes_GRch37_versions",
        "ensembl_genes_GRch37_location",
        "ensembl_genes_GRch37_ensembl_id",
        "ensembl_genes_GRch38_versions",
        "ensembl_genes_GRch38_location",
        "ensembl_genes_GRch38_ensembl_id",
        "hgnc_date_symbol_changed",
        "entity_type",
        "entity_name",
        "confidence_level",
        "penetrance",
        "mode_of_pathogenicity",
        "publications",
        "evidence",
        "phenotypes",
        "mode_of_inheritance",
        "tags",
        "panel_id",
        "panel_hash_id",
        "panel_name",
        "panel_disease_group",
        "panel_disease_sub_group",
        "panel_status",
        "panel_version",
        "panel_version_created",
        "panel_relevant_disorders",
        "panel_stats_number_of_genes",
        "panel_stats_number_of_strs",
        "panel_stats_number_of_regions",
        "panel_types_name",
        "panel_types_slug",
        "panel_types_description",
        "transcript",
        sep="\t",
        file=panelapp_output,
    )

    api_url = config["database"]["PANELAPP_AUS_API_URL"]
    # api_url = "https://panelapp.agha.umccr.org/api/v1/genes/"

    while api_url:
        response = requests.get(api_url)
        content = response.json()
        genes = content["results"]

        print(content["count"], api_url)

        for gene in genes:
            gene_data = gene["gene_data"]
            panel_info = gene["panel"]

            grch37_versions = []
            grch37_locations = []
            grch37_ensembl_ids = []
            if "GRch37" in gene_data["ensembl_genes"]:
                grch37_info = gene_data["ensembl_genes"]["GRch37"]
                for version, version_info in grch37_info.items():
                    grch37_versions.append(version)
                    grch37_locations.append(version_info.get("location", "-"))
                    grch37_ensembl_ids.append(
                        version_info.get("ensembl_id", "-")
                    )

            grch38_versions = []
            grch38_locations = []
            grch38_ensembl_ids = []
            if "GRch38" in gene_data["ensembl_genes"]:
                grch38_info = gene_data["ensembl_genes"]["GRch38"]
                for version, version_info in grch38_info.items():
                    grch38_versions.append(version)
                    grch38_locations.append(version_info.get("location", "-"))
                    grch38_ensembl_ids.append(
                        version_info.get("ensembl_id", "-")
                    )

            print(
                "||".join(gene_data["alias"]) if gene_data["alias"] else "-",
                gene_data["biotype"],
                gene_data["hgnc_id"],
                gene_data["gene_name"],
                "||".join(gene_data["omim_gene"])
                if gene_data["omim_gene"]
                else "-",
                "||".join(gene_data["alias_name"])
                if gene_data["alias_name"]
                else "-",
                gene_data["gene_symbol"],
                gene_data["hgnc_symbol"],
                gene_data["hgnc_release"],
                "||".join(grch37_versions) or "-",
                "||".join(grch37_locations) or "-",
                "||".join(grch37_ensembl_ids) or "-",
                "||".join(grch38_versions) or "-",
                "||".join(grch38_locations) or "-",
                "||".join(grch38_ensembl_ids) or "-",
                gene_data["hgnc_date_symbol_changed"],
                gene["entity_type"],
                gene["entity_name"],
                gene["confidence_level"],
                gene["penetrance"] or "-",
                gene["mode_of_pathogenicity"] or "-",
                "||".join(gene["publications"]).replace("\r\n", "")
                or "-",  # \r\n 있는 케이스가 1개 있었음
                "||".join(gene["evidence"]) or "-",
                "||".join(gene["phenotypes"]).replace("\t", " ") or "-",
                gene["mode_of_inheritance"],
                "||".join(gene["tags"]) or "-",
                panel_info["id"],
                panel_info["hash_id"],
                panel_info["name"],
                panel_info["disease_group"] or "-",
                panel_info["disease_sub_group"] or "-",
                panel_info["status"],
                panel_info["version"],
                panel_info["version_created"],
                "||".join(panel_info["relevant_disorders"]) or "-",
                panel_info["stats"]["number_of_genes"],
                panel_info["stats"]["number_of_strs"],
                panel_info["stats"]["number_of_regions"],
                "||".join(map(lambda x: x["name"], panel_info["types"])),
                "||".join(map(lambda x: x["slug"], panel_info["types"])),
                "||".join(map(lambda x: x["description"], panel_info["types"])),
                gene["transcript"] or "-",
                sep="\t",
                file=panelapp_output,
            )

        # api_url = None
        api_url = content["next"]

    panelapp_output.close()

    return


def parse_panelapp_json(panelapp_output_filepath, api_url):
    panelapp_output = open(panelapp_output_filepath, "w")

    print(
        "#alias",
        "biotype",
        "hgnc_id",
        "gene_name",
        "omim_gene",
        "alias_name",
        "gene_symbol",
        "hgnc_symbol",
        "hgnc_release",
        "ensembl_genes_GRch37_versions",
        "ensembl_genes_GRch37_location",
        "ensembl_genes_GRch37_ensembl_id",
        "ensembl_genes_GRch38_versions",
        "ensembl_genes_GRch38_location",
        "ensembl_genes_GRch38_ensembl_id",
        "hgnc_date_symbol_changed",
        "entity_type",
        "entity_name",
        "confidence_level",
        "penetrance",
        "mode_of_pathogenicity",
        "publications",
        "evidence",
        "phenotypes",
        "mode_of_inheritance",
        "tags",
        "panel_id",
        "panel_hash_id",
        "panel_name",
        "panel_disease_group",
        "panel_disease_sub_group",
        "panel_status",
        "panel_version",
        "panel_version_created",
        "panel_relevant_disorders",
        "panel_stats_number_of_genes",
        "panel_stats_number_of_strs",
        "panel_stats_number_of_regions",
        "panel_types_name",
        "panel_types_slug",
        "panel_types_description",
        "transcript",
        sep="\t",
        file=panelapp_output,
    )

    # api_url = "https://panelapp.genomicsengland.co.uk/api/v1/genes/?page=26"

    while api_url:
        response = requests.get(api_url)
        content = response.json()
        genes = content["results"]

        print(content["count"], api_url)

        for gene in genes:
            gene_data = gene["gene_data"]
            panel_info = gene["panel"]

            grch37_versions = []
            grch37_locations = []
            grch37_ensembl_ids = []
            if "GRch37" in gene_data["ensembl_genes"]:
                grch37_info = gene_data["ensembl_genes"]["GRch37"]
                for version, version_info in grch37_info.items():
                    grch37_versions.append(version)
                    grch37_locations.append(version_info.get("location", "-"))
                    grch37_ensembl_ids.append(
                        version_info.get("ensembl_id", "-")
                    )

            grch38_versions = []
            grch38_locations = []
            grch38_ensembl_ids = []
            if "GRch38" in gene_data["ensembl_genes"]:
                grch38_info = gene_data["ensembl_genes"]["GRch38"]
                for version, version_info in grch38_info.items():
                    grch38_versions.append(version)
                    grch38_locations.append(version_info.get("location", "-"))
                    grch38_ensembl_ids.append(
                        version_info.get("ensembl_id", "-")
                    )

            print(
                "||".join(gene_data["alias"]) if gene_data["alias"] else "-",
                gene_data["biotype"],
                gene_data["hgnc_id"],
                gene_data["gene_name"],
                "||".join(gene_data["omim_gene"])
                if gene_data["omim_gene"]
                else "-",
                "||".join(gene_data["alias_name"])
                if gene_data["alias_name"]
                else "-",
                gene_data["gene_symbol"],
                gene_data["hgnc_symbol"],
                gene_data["hgnc_release"],
                "||".join(grch37_versions) or "-",
                "||".join(grch37_locations) or "-",
                "||".join(grch37_ensembl_ids) or "-",
                "||".join(grch38_versions) or "-",
                "||".join(grch38_locations) or "-",
                "||".join(grch38_ensembl_ids) or "-",
                gene_data["hgnc_date_symbol_changed"],
                gene["entity_type"],
                gene["entity_name"],
                gene["confidence_level"],
                gene["penetrance"] or "-",
                gene["mode_of_pathogenicity"] or "-",
                "||".join(gene["publications"]).replace("\r\n", "")
                or "-",  # \r\n 있는 케이스가 1개 있었음
                "||".join(gene["evidence"]) or "-",
                "||".join(gene["phenotypes"]).replace("\t", " ") or "-",
                gene["mode_of_inheritance"],
                "||".join(gene["tags"]) or "-",
                panel_info["id"],
                panel_info["hash_id"],
                panel_info["name"],
                panel_info["disease_group"] or "-",
                panel_info["disease_sub_group"] or "-",
                panel_info["status"],
                panel_info["version"],
                panel_info["version_created"],
                "||".join(panel_info["relevant_disorders"]) or "-",
                panel_info["stats"]["number_of_genes"],
                panel_info["stats"]["number_of_strs"],
                panel_info["stats"]["number_of_regions"],
                "||".join(map(lambda x: x["name"], panel_info["types"])),
                "||".join(map(lambda x: x["slug"], panel_info["types"])),
                "||".join(map(lambda x: x["description"], panel_info["types"])),
                gene["transcript"] or "-",
                sep="\t",
                file=panelapp_output,
            )

        # api_url = None
        api_url = content["next"]

    panelapp_output.close()

    return


def parse_panelapp_rawdata(panelapp_rawdata_filepath: str, prefix: str) -> str:
    panelapp_parsed_filepath = f"panelapp.{prefix}.genes_parsed.tsv"
    panelapp_parsed_output = open(panelapp_parsed_filepath, "w")

    print(
        "#identifier",
        "hgnc_id",
        "hgnc_symbol",
        "confidence_level",
        "penetrance",
        "pmid",
        "omim_id",
        "inheritance",
        "panel_id",
        "panel_name",
        "panel_version",
        "tags",
        "line_number",
        sep="\t",
        file=panelapp_parsed_output,
    )

    ids = set()
    with open(panelapp_rawdata_filepath, "r") as panelapp:
        header = []
        line_number = 0
        for line in panelapp:
            line_number += 1
            if line.startswith("#"):
                header = line.strip().replace("#", "").split("\t")
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            hgnc_id = line_info["hgnc_id"]
            hgnc_symbol = line_info["hgnc_symbol"]
            confidence_level = line_info["confidence_level"]
            penetrance = line_info["penetrance"]
            publications = set(line_info["publications"].split("||"))
            phenotypes = line_info["phenotypes"].split("||")
            inheritance = line_info["mode_of_inheritance"]
            panel_id = line_info["panel_id"]
            panel_name = line_info["panel_name"]
            panel_version = line_info["panel_version"]
            tags = line_info["tags"]

            omim_ids = set()
            for phenotype in phenotypes:
                potential_omim_ids = set(
                    re.findall(r"\D\d\d\d\d\d\d\D", phenotype)
                )
                omim_ids |= potential_omim_ids

            if not omim_ids and publications == {"-"}:
                continue

            for omim_id in omim_ids:
                omim_id = omim_id[1:-1]
                identifier = hgnc_id.replace(":", "_") + "-" + "OMIM_" + omim_id

                ids.add(identifier)

                print(
                    identifier,
                    hgnc_id,
                    hgnc_symbol,
                    confidence_level,
                    penetrance,
                    "-",
                    omim_id,
                    inheritance,
                    panel_id,
                    panel_name,
                    panel_version,
                    tags,
                    line_number,
                    sep="\t",
                    file=panelapp_parsed_output,
                )

            if publications == {"-"}:
                continue

            for publication in publications:
                ####################################### 여기 publication preprocessing 필요함
                # re.sub(r"\(.*\)", "", publication)
                # pmid = re.match(r"\d+", publication)
                pmid = publication

                if not pmid:
                    continue

                identifier = hgnc_id.replace(":", "_") + "-" + "PMID_" + pmid

                ids.add(identifier)

                print(
                    identifier,
                    hgnc_id,
                    hgnc_symbol,
                    confidence_level,
                    penetrance,
                    pmid,
                    "-",
                    inheritance,
                    panel_id,
                    panel_name,
                    panel_version,
                    tags,
                    line_number,
                    sep="\t",
                    file=panelapp_parsed_output,
                )

    panelapp_parsed_output.close()
    print(len(ids))

    return panelapp_parsed_filepath


def remove_omim_from_panelapp(panelapp_file_path: str, prefix: str) -> str:
    omim_ids = get_omim_ids()

    panelapp_no_omim_filepath = f"panelapp.{prefix}.no_omim.tsv"
    panelapp_no_omim_output = open(panelapp_no_omim_filepath, "w")

    panelapp_omim_overlap_filepath = f"panelapp.{prefix}.omim_overlap.tsv"
    panelapp_omim_overlap_output = open(panelapp_omim_overlap_filepath, "w")

    with open(panelapp_file_path) as panelapp:
        header = []
        for line in panelapp:
            if line.startswith("#"):
                header = line.strip().replace("#", "").split("\t")
                print(line.strip(), file=panelapp_no_omim_output)
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            hgnc_omim_id = line_info["identifier"]

            if hgnc_omim_id not in omim_ids:
                print(line.strip(), file=panelapp_no_omim_output)
            else:
                print(line.strip(), file=panelapp_omim_overlap_output)

    panelapp_no_omim_output.close()
    panelapp_omim_overlap_output.close()

    return panelapp_no_omim_filepath


def remove_gencc_from_panelapp(
    panelapp_no_omim_filepath: str, prefix: str, keyword: str
) -> str:
    """_summary_

    Note:
        _description_

    Args:
        panelapp_filepath (str): _description_
        prefix (str): e.g. "uk", "aus"
        keyword (str): e.g. "Genomics England PanelApp", "PanelApp Australia"

    Returns:
        str: _description_

    Examples:
        >>> _description_
    """
    gencc_filepath = config["database"]["GENCC_PATH"]

    gencc_data = set()
    with open(gencc_filepath, "r") as gencc:
        header = []
        for line in gencc:
            if line.startswith("uuid"):
                header = line.strip().split("\t")
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            submitter = line_info["submitter_title"]
            hgnc_id = line_info["uuid"].split("-")[1]
            disease_id = line_info["uuid"].split("-")[2]

            if submitter == keyword:
                pmids = line_info["submitted_as_pmids"].split(", ")

                identifier = hgnc_id + "-" + disease_id
                gencc_data.add(identifier)

                for pmid in pmids:
                    if pmid == "-":
                        continue
                    identifier = hgnc_id + "-" + "PMID_" + pmid
                    gencc_data.add(identifier)

    panelapp_no_gencc_filepath = f"panelapp.{prefix}.no_omim.no_gencc.tsv"
    panelapp_no_gencc_output = open(panelapp_no_gencc_filepath, "w")

    with open(panelapp_no_omim_filepath, "r") as panelapp:
        header = []
        for line in panelapp:
            if line.startswith("#"):
                header = line.strip().replace("#", "").split("\t")
                print(line.strip(), file=panelapp_no_gencc_output)
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            identifier = line_info["identifier"]

            if identifier not in gencc_data:
                print(line.strip(), file=panelapp_no_gencc_output)

    panelapp_no_gencc_output.close()

    return panelapp_no_gencc_filepath


def aggregate_panelapp(panelapp_filepath, prefix):
    panelapp_info = {}
    with open(panelapp_filepath, "r") as panelapp:
        header = []
        for line in panelapp:
            if line.startswith("#"):
                header = line.strip().replace("#", "").split("\t")
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            hgnc_id = line_info["hgnc_id"]
            panel_id = line_info["panel_id"]
            panel_version = line_info["panel_version"]
            line_number = line_info["line_number"]

            new_id = (
                hgnc_id
                + "-"
                + panel_id
                + "-"
                + panel_version
                + "-"
                + line_number
            )

            if new_id in panelapp_info:
                panelapp_info[new_id]["pmid"].add(line_info["pmid"])
                panelapp_info[new_id]["omim_id"].add(line_info["omim_id"])
            else:
                panelapp_info[new_id] = {}
                panelapp_info[new_id]["hgnc_id"] = line_info["hgnc_id"]
                panelapp_info[new_id]["hgnc_symbol"] = line_info["hgnc_symbol"]
                panelapp_info[new_id]["confidence_level"] = line_info[
                    "confidence_level"
                ]
                panelapp_info[new_id]["penetrance"] = line_info["penetrance"]
                panelapp_info[new_id]["pmid"] = {line_info["pmid"]}
                panelapp_info[new_id]["omim_id"] = {line_info["omim_id"]}
                panelapp_info[new_id]["inheritance"] = line_info["inheritance"]
                panelapp_info[new_id]["panel_id"] = line_info["panel_id"]
                panelapp_info[new_id]["panel_name"] = line_info["panel_name"]
                panelapp_info[new_id]["panel_version"] = line_info[
                    "panel_version"
                ]
                panelapp_info[new_id]["line_number"] = line_info["line_number"]

    aggregate_panelapp_filepath = f"panelapp.{prefix}.aggregate.tsv"
    aggregate_panelapp_output = open(aggregate_panelapp_filepath, "w")
    print(
        "#hgnc_id",
        "hgnc_symbol",
        "confidence_level",
        "penetrance",
        "pmid",
        "omim_id",
        "inheritance",
        "panel_id",
        "panel_name",
        "panel_version",
        "line_number",
        sep="\t",
        file=aggregate_panelapp_output,
    )

    for content in panelapp_info.values():
        content["pmid"].discard("-")
        content["omim_id"].discard("-")

        print(
            content["hgnc_id"],
            content["hgnc_symbol"],
            content["confidence_level"],
            content["penetrance"],
            "||".join(content["pmid"]) or "-",
            "||".join(content["omim_id"]) or "-",
            content["inheritance"],
            content["panel_id"],
            content["panel_name"],
            content["panel_version"],
            content["line_number"],
            sep="\t",
            file=aggregate_panelapp_output,
        )

    aggregate_panelapp_output.close()

    return


def get_highest_panels(panelapp_filepath: str, prefix: str) -> str:
    panelapp_filtered_filepath = (
        f"panelapp.{prefix}.aggregate.highest_panels.tsv"
    )
    panelapp_filtered_output = open(panelapp_filtered_filepath, "w")

    panelapp_filtered_data = {}
    with open(panelapp_filepath, "r") as panelapp:
        header = []
        for line in panelapp:
            if line.startswith("#"):
                header = line.strip().replace("#", "").split("\t")
                print(line.strip(), file=panelapp_filtered_output)
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            hgnc_id = line_info["hgnc_id"]
            panel_id = line_info["panel_id"]
            panel_major_version, panel_minor_version = map(
                int, line_info["panel_version"].split(".")
            )

            if panel_major_version < 1:
                continue

            key = hgnc_id + "-" + panel_id

            if key in panelapp_filtered_data:
                if panel_major_version > panelapp_filtered_data[key][
                    "panel_major_version"
                ] or (
                    panel_major_version
                    == panelapp_filtered_data[key]["panel_major_version"]
                    and panel_minor_version
                    > panelapp_filtered_data[key]["panel_minor_version"]
                ):
                    panelapp_filtered_data[key][
                        "panel_major_version"
                    ] = panel_major_version
                    panelapp_filtered_data[key][
                        "panel_minor_version"
                    ] = panel_minor_version
                    panelapp_filtered_data[key]["content"] = line.strip()
            else:
                panelapp_filtered_data[key] = {}
                panelapp_filtered_data[key][
                    "panel_major_version"
                ] = panel_major_version
                panelapp_filtered_data[key][
                    "panel_minor_version"
                ] = panel_minor_version
                panelapp_filtered_data[key]["content"] = line.strip()

    for hgnc_id, filtered_info in panelapp_filtered_data.items():
        print(filtered_info["content"], file=panelapp_filtered_output)

    panelapp_filtered_output.close()

    return panelapp_filtered_filepath


def filter_panelapp_aggregate(panelapp_filepath: str, prefix: str) -> str:
    panelapp_filtered_filepath = f"panelapp.{prefix}.aggregate.filtered.tsv"
    # panelapp_filtered_filepath = (
    #     f"panelapp.{prefix}.aggregate.filtered_only_omim.tsv"
    # )
    panelapp_filtered_output = open(panelapp_filtered_filepath, "w")

    with open(panelapp_filepath, "r") as panelapp:
        header = []
        for line in panelapp:
            if line.startswith("#"):
                header = line.strip().replace("#", "").split("\t")
                print(line.strip(), file=panelapp_filtered_output)
                continue

            columns = line.strip().split("\t")
            line_info = dict(zip(header, columns))

            pmid = line_info["pmid"]
            omim_id = line_info["omim_id"]

            pmids = pmid.split("||")

            if omim_id != "-":
                print(line.strip(), file=panelapp_filtered_output)
            else:
                if len(pmids) >= 2:
                    print(line.strip(), file=panelapp_filtered_output)

    panelapp_filtered_output.close()

    return panelapp_filtered_filepath


def main():
    # remove_omim_from_gencc()

    panelapp_uk_filepath = config["database"]["PANELAPP_UK_PATH"]
    panelapp_uk_api = config["database"]["PANELAPP_UK_API_URL"]
    panelapp_aus_filepath = config["database"]["PANELAPP_AUS_PATH"]
    panelapp_aus_api = config["database"]["PANELAPP_AUS_API_URL"]

    # parse_panelapp_json(panelapp_uk_filepath, panelapp_uk_api)
    # parse_panelapp_rawdata(panelapp_uk_filepath, "uk")
    # remove_omim_from_panelapp("panelapp.uk.genes_parsed.tsv", "uk")
    # remove_gencc_from_panelapp(
    #     "panelapp.uk.no_omim.tsv", "uk", "Genomics England PanelApp"
    # )
    # aggregate_panelapp("panelapp.uk.no_omim.no_gencc.tsv", "uk")
    # get_highest_panels("panelapp.uk.aggregate.tsv", "uk")
    filter_panelapp_aggregate("panelapp.uk.aggregate.highest_panels.tsv", "uk")

    # parse_panelapp_json(panelapp_aus_filepath, panelapp_aus_api)
    # parse_panelapp_rawdata(panelapp_aus_filepath, "aus")
    # remove_omim_from_panelapp("panelapp.aus.genes_parsed.tsv", "aus")
    # remove_gencc_from_panelapp(
    #     "panelapp.aus.no_omim.tsv", "aus", "PanelApp Australia"
    # )
    # aggregate_panelapp("panelapp.aus.no_omim.no_gencc.tsv", "aus")
    # get_highest_panels("panelapp.aus.aggregate.tsv", "aus")
    filter_panelapp_aggregate(
        "panelapp.aus.aggregate.highest_panels.tsv", "aus"
    )

    # parse_panelapp_aus_json()


if __name__ == "__main__":
    main()
