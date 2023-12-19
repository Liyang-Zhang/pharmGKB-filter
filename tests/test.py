import json

import pandas as pd

config_path = (
    "/home/leon/repo/pharmGKB-extraction/pharmgkb_extraction/config/custom_class.json"
)

with open(config_path, "r") as config_file:
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

clinical_annotations_path = "/home/leon/data/lung65/chemo/clinical_annotations.tsv"

df_cli_anno = pd.read_csv(clinical_annotations_path, sep="\t", header=0)

# Apply all conditions in a single line
# filtered_df = df_cli_anno[
#    df_cli_anno["Drug(s)"].apply(
#        lambda x: isinstance(x, str) and any(drug in x.lower() for drug in drug_list)
#    )
#    & df_cli_anno["Gene"].apply(
#        lambda x: isinstance(x, str) and any(gene in x.lower() for gene in gene_list)
#    )
#    &
#    df_cli_anno["Phenotype(s)"].apply(
#        lambda x: isinstance(x, str)
#        and any(phenotype in x for phenotype in phenotype_list)
#    )
#    & df_cli_anno["Level of Evidence"].apply(
#        lambda x: isinstance(x, str) and any(level in x for level in level_list)
#    )
# ]


def any_in_string(items, string):
    return any(item in string.lower() for item in items)


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


# Output the filtered DataFrame to a new TSV file
output_path = "/home/leon/repo/pharmGKB-extraction/tests/lung65_clinical_annotation.tsv"

df_cli_anno.to_csv(output_path, sep="\t", index=False)

print(f"Filtered data has been saved to {output_path}")
