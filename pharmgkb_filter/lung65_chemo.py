import json
import logging
from pathlib import Path

import pandas as pd
from intervaltree import IntervalTree
from pysam import VariantFile

from pharmgkb_filter.logging.helper import LoggingLevel, config_console_debug_logger

log = logging.getLogger(__name__)
config_console_debug_logger(__name__)

logging.basicConfig(level=LoggingLevel.INFO.upper())


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

    log.info(f"Gene list:{gene_list}")
    log.info(f"Drug list: {drug_list}")
    log.info(f"Phenotype list: {phenotype_list}")
    log.info(f"Level list: {level_list}")

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


def add_hg19_genome_location_info(df: pd.DataFrame, vcf: Path) -> pd.DataFrame:
    """
    Input:
        1. pharmGKB clinical annotations dataframe
        2. matched result VCF file
    """
    # generate a genome-to-rs dataframe for merging
    pos_to_rs_list = []
    with VariantFile(str(vcf), "r") as vcf_in:
        for record in vcf_in:
            pos_to_rs_list.append(
                {"chrom": record.chrom, "pos": record.pos, "id": record.id}
            )
    pos_to_rs_df = pd.DataFrame(pos_to_rs_list)

    merged_df = df.merge(
        pos_to_rs_df, left_on="Variant/Haplotypes", right_on="id", how="left"
    )

    merged_df[["chrom", "pos", "id"]] = merged_df[["chrom", "pos", "id"]].fillna(".")

    return merged_df


def add_interval_column(
    df: pd.DataFrame, bed_file: Path, column_name: str
) -> pd.DataFrame:
    # Read the BED file and create an interval tree for each chromosome
    interval_trees = {}
    with open(bed_file, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom, start, end = parts[:3]
                start, end = int(start), int(end)
                if chrom not in interval_trees:
                    interval_trees[chrom] = IntervalTree()
                interval_trees[chrom][start:end] = True

    # Function to check if a position is in any interval
    def is_in_interval(row):
        chrom = row["chrom"]
        pos = row["pos"]
        if chrom in interval_trees and interval_trees[chrom][pos]:
            return True
        return False

    # Apply the function to each row
    df[column_name] = df.apply(is_in_interval, axis=1)

    return df


def main():
    # fill the path before running
    config_path = Path()
    clinical_annotations_path = Path()
    rs_matched_vcf_path = Path()
    bed_133 = Path()
    bed_188 = Path()

    output_path = ()

    df_cli_anno = pd.read_csv(clinical_annotations_path, sep="\t", header=0)
    df_cli_anno_common_filters = add_common_filter_columns(df_cli_anno, config_path)

    df_cli_anno_pos_merge = add_hg19_genome_location_info(
        df_cli_anno_common_filters, rs_matched_vcf_path
    )

    df_cli_anno_pos_bed_check = add_interval_column(
        df_cli_anno_pos_merge, bed_133, "is_in_133_interval"
    )

    df_cli_anno_pos_bed_check = add_interval_column(
        df_cli_anno_pos_bed_check, bed_188, "is_in_188_interval"
    )

    df_cli_anno_pos_bed_check.to_csv(output_path, sep="\t", index=False)
    print(f"Filtered data has been saved to {output_path}")


if __name__ == "__main__":
    main()
