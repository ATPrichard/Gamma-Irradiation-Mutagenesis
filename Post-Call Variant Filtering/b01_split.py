#!/usr/bin/env python3
import os
import sys


# CONFIG SECTION: 

# INPUT
INPUT_VCF = "vcf/GB2_call.vcf"

# OUTPUT
SNP_VCF = "Results/snps/GB2_raw.snps.vcf"
INDEL_VCF = "Results/indels/GB2_raw.indels.vcf"

# SUMMARY PATH
SUMMARY_FILE = "Summaries/GB2_split_summary.txt"

# SNP
SNP_REF_LENGTH = 1
SNP_ALT_LENGTH = 1

SUMMARY_PRECISION = 4

# Progress Tracking
PROGRESS_INTERVAL = 50000


def is_snp(ref, alts):
    if len(ref) != SNP_REF_LENGTH:
        return False
    for alt in alts:
        if len(alt) != SNP_ALT_LENGTH:
            return False
    return True


def ensure_dirs(*paths):
    for path in paths:
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)


def count_variants(vcf_path):
    """
    Count non-header lines for progress tracking.
    """
    count = 0
    with open(vcf_path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    return count


def split_vcf():
    try:
        total_variants = count_variants(INPUT_VCF)
    except FileNotFoundError:
        sys.exit(f"ERROR: Input VCF not found: {INPUT_VCF}")

    processed = 0
    snp_count = 0
    indel_count = 0

    with open(INPUT_VCF, "r") as infile, \
         open(SNP_VCF, "w") as snp_out, \
         open(INDEL_VCF, "w") as indel_out:

        for line in infile:
            if line.startswith("#"):
                snp_out.write(line)
                indel_out.write(line)
                continue

            processed += 1
            fields = line.rstrip().split("\t")
            ref = fields[3]
            alts = fields[4].split(",")

            if is_snp(ref, alts):
                snp_out.write(line)
                snp_count += 1
            else:
                indel_out.write(line)
                indel_count += 1

            if processed % PROGRESS_INTERVAL == 0 or processed == total_variants:
                percent = (processed / total_variants) * 100
                print(
                    f"[Step 0] Processed {processed:,}/{total_variants:,} "
                    f"variants ({percent:.1f}%)",
                    flush=True
                )

    with open(SUMMARY_FILE, "w") as s:
        s.write("Step 0: Split SNPs and INDELs\n")
        s.write("============================\n\n")
        s.write(f"Input VCF: {INPUT_VCF}\n\n")
        s.write(f"Total variants: {processed}\n")
        s.write(f"SNPs: {snp_count}\n")
        s.write(f"INDELs: {indel_count}\n")

        if processed > 0:
            s.write(
                f"SNP fraction: {snp_count / processed:.{SUMMARY_PRECISION}f}\n"
            )
            s.write(
                f"INDEL fraction: {indel_count / processed:.{SUMMARY_PRECISION}f}\n"
            )


def main():
    ensure_dirs(SNP_VCF, INDEL_VCF, SUMMARY_FILE)
    split_vcf()


if __name__ == "__main__":
    main()
