#!/usr/bin/env python3

import os
import re
import sys

#CONFIG SECTION:

SNP_VCF = "SampleX_SNP.vcf"
INDEL_VCF = "SampleX_INDEL.vcf"
OUTPUT_VCF = "SampleX_FINAL.vcf"

STRICT_HEADER_CHECK = True


def read_vcf(path):
    meta_header = []
    column_header = None
    variants = []

    with open(path, "r") as f:
        for line in f:
            if line.startswith("##"):
                meta_header.append(line)
            elif line.startswith("#CHROM"):
                column_header = line
            else:
                if line.strip():
                    variants.append(line)

    if column_header is None:
        raise ValueError(f"No #CHROM header line found in {path}")

    return meta_header, column_header, variants


def chrom_sort_key(chrom):
    c = chrom.strip()

    match = re.match(r"NC_(\d+)\.\d+$", c)
    if match:
        return (0, int(match.group(1)))

    if c.lower().startswith("chr"):
        c = c[3:]

    if c.isdigit():
        return (1, int(c))

    c_upper = c.upper()
    if c_upper == "X":
        return (2, 23)
    if c_upper == "Y":
        return (2, 24)
    if c_upper in ("M", "MT"):
        return (2, 25)

    return (3, c_upper)


def variant_sort_key(line):
    fields = line.rstrip("\n").split("\t")

    if len(fields) < 2:
        raise ValueError(f"Malformed variant line: {line.strip()}")

    chrom = fields[0]

    try:
        pos = int(fields[1])
    except ValueError:
        raise ValueError(f"Invalid POS value in line: {line.strip()}")

    return (chrom_sort_key(chrom), pos)


def main():
    if not os.path.exists(SNP_VCF):
        print(f"ERROR: SNP VCF not found: {SNP_VCF}")
        sys.exit(1)

    if not os.path.exists(INDEL_VCF):
        print(f"ERROR: Indel VCF not found: {INDEL_VCF}")
        sys.exit(1)

    print("Reading SNP VCF...")
    snp_meta, snp_header, snp_variants = read_vcf(SNP_VCF)

    print("Reading Indel VCF...")
    indel_meta, indel_header, indel_variants = read_vcf(INDEL_VCF)


    if snp_header != indel_header:
        print("ERROR: #CHROM header lines do not match between files.")
        print("SNP header:  ", snp_header.strip())
        print("Indel header:", indel_header.strip())
        sys.exit(1)


    if STRICT_HEADER_CHECK and snp_meta != indel_meta:
        print("ERROR: Metadata header lines (##) do not match between files.")
        print("Set STRICT_HEADER_CHECK = False if you want to ignore this.")
        sys.exit(1)

    print(f"SNP variants:   {len(snp_variants)}")
    print(f"Indel variants: {len(indel_variants)}")

    combined_variants = snp_variants + indel_variants

    print("Sorting combined variants by chromosome and position...")
    combined_variants.sort(key=variant_sort_key)

    print(f"Writing output: {OUTPUT_VCF}")
    with open(OUTPUT_VCF, "w") as out:
        for line in snp_meta:
            out.write(line)
        out.write(snp_header)
        for line in combined_variants:
            out.write(line)

    print("Done.")
    print(f"Combined variant count: {len(combined_variants)}")


if __name__ == "__main__":
    main()
