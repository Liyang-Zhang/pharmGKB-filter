import json
from pathlib import Path

import pandas as pd


def any_in_string(items, string):
    return any(item in string.lower() for item in items)


def add_common_filter_columns(
    df_cli_anno: pd.DataFrame, config_path: Path
) -> pd.DataFrame:
    with open(str(config_path), "r") as config_file:
        data = json.load(config_file)

    gene_list = [gene.lower() for gene in data["gene"]]
    drug_list = [drug.lower() for drug in data["drug"]]
    phenotype_list_ld = [drug.lower() for drug in data["phenotype_ld"]]
    phenotype_list = [phenotype.lower() for phenotype in data["phenotype"]]
    level_list = data["level"]

    print(gene_list)
    print(drug_list)
    print(phenotype_list)
    print(level_list)

    df_cli_anno["phenotype_ld_1219"] = df_cli_anno["Phenotype(s)"].apply(
        lambda x: any_in_string(phenotype_list_ld, x) if isinstance(x, str) else False
    )
    df_cli_anno["Matches Drug"] = df_cli_anno["Drug(s)"].apply(
        lambda x: any_in_string(drug_list, x) if isinstance(x, str) else False
    )
    df_cli_anno["Matches Gene"] = df_cli_anno["Gene"].apply(
        lambda x: any_in_string(gene_list, x) if isinstance(x, str) else False
    )
    df_cli_anno["Matches Phenotype"] = df_cli_anno["Phenotype(s)"].apply(
        lambda x: any_in_string(phenotype_list, x) if isinstance(x, str) else False
    )
    df_cli_anno["Matches Level"] = df_cli_anno["Level of Evidence"].apply(
        lambda x: x in level_list if isinstance(x, str) else False
    )

    return df_cli_anno


# TODO Add genome location to each row
# Input: 1. matched result VCF file or bed file
def add_hg19_genome_location_info():
    ...


# TODO Add columns to tell if the record rs ID is in the designed probe intervals
# Input: 1. porbe ned files
def add_if_covered_by_probe_info():
    ...


def main():
    config_path = Path(
        "/Users/zhangliyang/repo/pharmGKB-extraction/pharmgkb_extraction/config/custom_class.json"
    )
    clinical_annotations_path = Path(
        "/Users/zhangliyang/repo/pharmGKB-extraction/tests/clinicalAnnotations/clinical_annotations.tsv"
    )

    df_cli_anno = pd.read_csv(clinical_annotations_path, sep="\t", header=0)
    df_cli_anno = add_common_filter_columns(df_cli_anno, config_path)

    output_path = "/Users/zhangliyang/repo/pharmGKB-extraction/tests/clinicalAnnotations/lung65_clinical_annotation.tsv"
    df_cli_anno.to_csv(output_path, sep="\t", index=False)

    print(f"Filtered data has been saved to {output_path}")


if __name__ == "__main__":
    main()
