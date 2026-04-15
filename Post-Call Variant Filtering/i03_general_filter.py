#!/usr/bin/env python3

import os
import sys

# CONFIG SECTION

INPUT_LIST = "VCF_lists/GB2_GF_list.txt"
RESULTS_DIR = "Results"
SUMMARY_DIR = "Summaries"

MIN_IDV = 15
MIN_IMF = 0.2
MIN_STRAND_READS = 3  

PROGRESS_UPDATE_INTERVAL = 25

def parse_info_field(info_field):
    info_dict = {}
    for entry in info_field.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_dict[key] = value
    return info_dict

def safe_basename_noext(path):
    base = os.path.basename(path)
    if base.lower().endswith(".vcf"):
        return base[:-4]
    return os.path.splitext(base)[0]

def is_heterozygous(gt):
    return gt in {"0/1", "1/0", "0|1", "1|0"}


def process_single_vcf(vcf_path):
    if not os.path.isfile(vcf_path):
        print(f"WARNING: input VCF not found, skipping: {vcf_path}", file=sys.stderr)
        return

    input_base = safe_basename_noext(vcf_path)
    output_vcf = os.path.join(RESULTS_DIR, f"{input_base}_GF.vcf")
    summary_file = os.path.join(SUMMARY_DIR, f"{input_base}_GF_summary.txt")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(SUMMARY_DIR, exist_ok=True)

    kept = 0
    discarded = 0
    total = 0
    strand_failed = 0

    with open(vcf_path, "r") as infile, open(output_vcf, "w") as out:
        headers = []

        # Header processing
        for line in infile:
            if line.startswith("##"):
                headers.append(line)
            elif line.startswith("#CHROM"):
                headers.append(line)
                out.writelines(headers)
                header_parts = line.strip().split("\t")
                sample_columns = header_parts[9:]
                break

        # Variants
        for i, line in enumerate(infile, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 10:
                discarded += 1
                total += 1
                continue

            info_dict = parse_info_field(parts[7])
            try:
                idv = int(info_dict.get("IDV", 0))
                imf = float(info_dict.get("IMF", 0))
            except ValueError:
                discarded += 1
                total += 1
                continue

            if idv < MIN_IDV or imf < MIN_IMF:
                discarded += 1
                total += 1
                continue

            # FORMAT parsing for all samples
            fmt_keys = parts[8].split(":")
            all_pass = True

            for sample_val in parts[9:]:
                sample_data = sample_val.split(":")
                if len(sample_data) != len(fmt_keys):
                    all_pass = False
                    strand_failed += 1
                    break
                fmt = dict(zip(fmt_keys, sample_data))
                gt = fmt.get("GT", "")

                if is_heterozygous(gt):
                    try:
                        adf = list(map(int, fmt["ADF"].split(",")))
                        adr = list(map(int, fmt["ADR"].split(",")))

                        # Must have 2 alleles per strand
                        if len(adf) < 2 or len(adr) < 2:
                            all_pass = False
                            strand_failed += 1
                            break

                        # Strict filter
                        if any(x < MIN_STRAND_READS for x in (adf[0], adf[1], adr[0], adr[1])):
                            all_pass = False
                            strand_failed += 1
                            break

                    except (ValueError, KeyError):
                        all_pass = False
                        strand_failed += 1
                        break

            if all_pass:
                out.write(line)
                kept += 1
            else:
                discarded += 1

            total += 1
            if i % PROGRESS_UPDATE_INTERVAL == 0:
                print(f"[{input_base}] Processed {i:,} lines...")

    percent_kept = (kept / total * 100) if total else 0
    percent_discarded = (discarded / total * 100) if total else 0

    with open(summary_file, "w") as s:
        s.write(f"Input VCF: {vcf_path}\n")
        s.write(f"Output VCF: {output_vcf}\n\n")

        s.write("Filter parameters:\n")
        s.write(f"  MIN_IDV = {MIN_IDV}\n")
        s.write(f"  MIN_IMF = {MIN_IMF}\n")
        s.write(
            f"  Heterozygous strand support: "
            f">= {MIN_STRAND_READS} reads per allele per strand (ADF/ADR)\n"
        )
        s.write(f"  Applied to all samples in the VCF\n\n")

        s.write(f"Total sites processed: {total}\n")
        s.write(f"Sites kept: {kept}\n")
        s.write(f"Sites discarded: {discarded}\n")
        s.write(f"  - Failed strand support: {strand_failed}\n\n")
        s.write(f"Percent kept: {percent_kept:.2f}%\n")
        s.write(f"Percent discarded: {percent_discarded:.2f}%\n")

    print(f"Finished {vcf_path}: kept {kept}, discarded {discarded}")
    print(f"Summary -> {summary_file}\n")


def main():
    if not os.path.isfile(INPUT_LIST):
        print(f"ERROR: INPUT_LIST not found: {INPUT_LIST}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(SUMMARY_DIR, exist_ok=True)

    with open(INPUT_LIST, "r") as f:
        vcf_paths = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]

    for vcf_path in vcf_paths:
        process_single_vcf(vcf_path)

    print("All VCFs processed successfully!")


if __name__ == "__main__":
    main()

