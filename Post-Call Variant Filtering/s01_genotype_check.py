#!/usr/bin/env python3

import os

# CONFIG SECTION
INPUT_VCF = "Results/snps/GB2_raw.snps.vcf"

RESULTS_DIR = "Results"
SUMMARY_DIR = "Summaries"

PROGRESS_UPDATE_INTERVAL = 100000

MISSING_THRESHOLD = 70

VALID_GTS = {"0/0", "0/1", "1/1", "./."}


def parse_genotype(field):
    parts = field.split(":")
    gt = parts[0]
    ad = parts[4] if len(parts) > 4 else None
    return gt, ad, parts

def check_homozygous_purity(gt, ad_field):
    if gt not in ("0/0", "1/1"):
        return True

    if ad_field is None or ad_field == ".":
        return True 

    try:
        ref, alt = map(int, ad_field.split(",")[:2])
    except ValueError:
        return True

    if gt == "0/0":
        return alt == 0
    if gt == "1/1":
        return ref == 0

    return True

def process_vcf(input_vcf):

    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(SUMMARY_DIR, exist_ok=True)

    base_name = os.path.basename(input_vcf).split(".")[0]
    output_vcf = os.path.join(RESULTS_DIR, f"{base_name}_GT.vcf")
    summary_file = os.path.join(SUMMARY_DIR, f"{base_name}_GT_summary.txt")

    total_sites = 0
    sites_removed = 0
    sites_with_missing = 0
    total_masked = 0

    with open(input_vcf, "r") as infile, open(output_vcf, "w") as out:

        headers = []
        sample_names = []
        per_sample_masked = []

        # Read headers
        for line in infile:
            if line.startswith("##"):
                headers.append(line)
            elif line.startswith("#CHROM"):
                headers.append(line)
                header_parts = line.strip().split("\t")
                sample_names = header_parts[9:]
                per_sample_masked = [0] * len(sample_names)
                out.writelines(headers)
                break

        # Process variant lines
        for i, line in enumerate(infile, start=1):

            if not line.strip() or line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            alt_field = parts[4]

            total_sites += 1

            # Remove multi-allelic sites
            if "," in alt_field:
                sites_removed += 1
                continue

            sample_fields = parts[9:]
            new_sample_fields = []

            remove_site = False
            missing_count = 0

            for idx, field in enumerate(sample_fields):

                gt, ad, field_parts = parse_genotype(field)

                # Remove site if GT not allowed
                if gt not in VALID_GTS:
                    remove_site = True
                    break

                # Mask failing homozygous purity
                if not check_homozygous_purity(gt, ad):
                    field_parts[0] = "./."
                    total_masked += 1
                    per_sample_masked[idx] += 1
                    gt = "./."

                if gt == "./.":
                    missing_count += 1

                new_sample_fields.append(":".join(field_parts))

            if remove_site:
                sites_removed += 1
                continue

            missing_pct = (missing_count / len(sample_fields)) * 100

            if missing_count >= 1:
                sites_with_missing += 1

            if missing_pct >= MISSING_THRESHOLD:
                sites_removed += 1
                continue

            parts[9:] = new_sample_fields
            out.write("\t".join(parts) + "\n")

            if i % PROGRESS_UPDATE_INTERVAL == 0:
                print(f"Processed {i:,} lines...")


    with open(summary_file, "w") as s:
        s.write(f"Total sites processed: {total_sites}\n")
        s.write(f"Genotypes masked: {total_masked}\n")
        s.write(f"Sites removed: {sites_removed}\n")
        s.write(f"Sites with >=1 ./.: {sites_with_missing}\n\n")
        s.write("Per-sample masked counts:\n")
        for name, count in zip(sample_names, per_sample_masked):
            s.write(f"  {name}: {count}\n")

    print("Finished processing.")
    print(f"Summary written to: {summary_file}")


if __name__ == "__main__":
    if not os.path.isfile(INPUT_VCF):
        print(f"Input VCF not found: {INPUT_VCF}")
    else:
        process_vcf(INPUT_VCF)
