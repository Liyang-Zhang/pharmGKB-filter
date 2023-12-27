import gzip
import json
from pathlib import Path


# Convert vcf to dicts by stream processing
def process_vcf_stream(gzipped_vcf_file_path: Path, gzipped_json_file_path: Path):
    """
    Converts VCF file to JSON format.
    """
    with gzip.open(gzipped_vcf_file_path, "rt") as vcf_file, gzip.open(
        gzipped_json_file_path, "wt"
    ) as output_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) > 2:
                key = fields[2]
                value = line.strip()
                output_file.write(json.dumps({key: value}) + "\n")


# Convert all dicts into single dict with multiple json file
def batch_process_json_objects(
    gzipped_json_file_path: Path, batch_size: int, output_file_path: Path
):
    """
    Processes JSON objects in batches.
    """
    batch_count = 0
    current_batch = {}

    with gzip.open(gzipped_json_file_path, "rt") as file:
        for line in file:
            json_obj = json.loads(line)
            current_batch.update(json_obj)

            if len(current_batch) >= batch_size:
                write_batch_to_file(
                    current_batch, f"{output_file_path}_{batch_count}.json.gz"
                )
                current_batch = {}
                batch_count += 1

        if current_batch:
            write_batch_to_file(
                current_batch, f"{output_file_path}_{batch_count}.json.gz"
            )


def write_batch_to_file(batch_data: dict, file_name: str):
    with gzip.open(file_name, "wt") as file:
        json.dump(batch_data, file)


# Extract vcf record in a given rs list in all sub json files
def process_all_sub_json_files(json_dir: Path, rs_list: Path, output_vcf: Path):
    """
    Processes sub JSON files to extract specific VCF records.
    """
    for json_file in json_dir.glob(
        "*.json.gz"
    ):  # Iterate over all .json.gz files in the directory
        print(json_file)
        with gzip.open(json_file, "rt") as file:
            vcf_dict = json.load(file)

        with open(rs_list, "r") as rs_file, open(
            output_vcf, "w"
        ) as out_f:  # Append mode for output file
            for rs in rs_file:
                rs = rs.strip()  # Strip any whitespace or newline characters
                if rs in vcf_dict:
                    print(f"{rs} found in {json_file}")
                    rs_record = vcf_dict[rs]
                    out_f.write(rs_record + "\n")


gzipped_vcf_file_path = Path(
    "/home/leon/data/lung65/dbSNP_latest_release/GCF_000001405.25.gz"
)

gzipped_json_file_path = Path(
    "/home/leon/data/lung65/dbSNP_latest_release/GCF_000001405.25.json.gz"
)

batch_dir = Path("/home/leon/data/lung65/dbSNP_latest_release/GCF_000001405.25.merged")

rs_list = Path(
    "/home/leon/data/lung65/dbSNP_latest_release/uniq_VariantHaplotypes_list.txt"
)

json_dir = Path("/home/leon/data/lung65/dbSNP_latest_release/dbSNP156_sub_jsons")

output_vcf = Path("/home/leon/data/lung65/dbSNP_latest_release/test.out.vcf")

# Time consuming, should run one function a time
process_vcf_stream(gzipped_vcf_file_path, gzipped_json_file_path)
batch_process_json_objects(gzipped_json_file_path, 10000000, batch_dir)
process_all_sub_json_files(json_dir, rs_list, output_vcf)
