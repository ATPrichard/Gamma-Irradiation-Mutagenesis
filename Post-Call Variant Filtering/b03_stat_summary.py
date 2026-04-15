#!/usr/bin/env python3

import os
from collections import defaultdict

# CONFIG Section

# List of VCF final VCFs
INPUT_VCF_LIST_FILE = "Final_VCF_list.txt"

OUTPUT_DIR = "Final VCF/"
COMBINED_FILE = "final_summary.txt"

GENOME_SIZE_BP = 133_198_052
CHROM_SIZE_FILE = "Text Files/KAP4_Chrom_Sizes.txt"


def load_vcf_paths():

    paths = list(INPUT_VCFS)

    if INPUT_VCF_LIST_FILE and os.path.exists(INPUT_VCF_LIST_FILE):
        with open(INPUT_VCF_LIST_FILE, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                paths.append(line)


    seen = set()
    unique_paths = []
    for p in paths:
        if p not in seen:
            unique_paths.append(p)
            seen.add(p)

    return unique_paths


def load_chrom_sizes(chrom_size_file):
    chrom_sizes = {}
    if chrom_size_file and os.path.exists(chrom_size_file):
        with open(chrom_size_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                chrom = parts[0]
                size = int(parts[1])
                chrom_sizes[chrom] = size
    return chrom_sizes


def is_snp(ref, alt_alleles):
    return len(ref) == 1 and all(len(a) == 1 for a in alt_alleles)


def is_indel(ref, alt_alleles):
    return any(len(a) != len(ref) for a in alt_alleles)


def is_transition(ref, alt):
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    return (ref, alt) in transitions


def normalize_substitution(ref, alt):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    ref = ref.upper()
    alt = alt.upper()

    if ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"}:
        return None

    if ref in {"C", "T"}:
        return f"{ref}>{alt}"

    ref_rc = complement[ref]
    alt_rc = complement[alt]
    return f"{ref_rc}>{alt_rc}"


def parse_vcf(vcf_path, chrom_sizes=None, genome_size_bp=None):
    total_variants = 0
    snp_count = 0
    indel_count = 0
    qual_sum = 0.0
    indel_sizes = []
    chrom_counts = defaultdict(lambda: {"snps": 0, "indels": 0, "total": 0})
    ti = tv = 0

    insertion_count = 0
    deletion_count = 0
    multiallelic_count = 0

    mutation_spectrum = {
        "C>A": 0,
        "C>G": 0,
        "C>T": 0,
        "T>A": 0,
        "T>C": 0,
        "T>G": 0,
    }

    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue

            chrom, _, _, ref, alt, qual = parts[:6]
            alt_alleles = alt.split(",")
            qual_val = float(qual) if qual != "." else 0.0

            total_variants += 1
            qual_sum += qual_val
            chrom_counts[chrom]["total"] += 1

            if len(alt_alleles) > 1:
                multiallelic_count += 1

            if is_snp(ref, alt_alleles):
                snp_count += 1
                chrom_counts[chrom]["snps"] += 1

                for a in alt_alleles:
                    if is_transition(ref, a):
                        ti += 1
                    else:
                        tv += 1

                    mut_class = normalize_substitution(ref, a)
                    if mut_class in mutation_spectrum:
                        mutation_spectrum[mut_class] += 1

            elif is_indel(ref, alt_alleles):
                indel_count += 1
                chrom_counts[chrom]["indels"] += 1

                for a in alt_alleles:
                    size_diff = len(a) - len(ref)
                    indel_sizes.append(abs(size_diff))

                    if size_diff > 0:
                        insertion_count += 1
                    elif size_diff < 0:
                        deletion_count += 1

    avg_qual = qual_sum / total_variants if total_variants else 0
    avg_snps_per_chr = snp_count / len(chrom_counts) if chrom_counts else 0
    avg_indels_per_chr = indel_count / len(chrom_counts) if chrom_counts else 0
    largest_indel = max(indel_sizes) if indel_sizes else 0
    mean_indel_size = (sum(indel_sizes) / len(indel_sizes)) if indel_sizes else 0
    titv_ratio = round(ti / tv, 3) if tv > 0 else "NA"

    snp_prop = round((snp_count / total_variants) * 100, 2) if total_variants else 0
    indel_prop = round((indel_count / total_variants) * 100, 2) if total_variants else 0
    snp_indel_ratio = round(snp_count / indel_count, 3) if indel_count > 0 else "NA"

    # Overall variant density
    if genome_size_bp and genome_size_bp > 0:
        variant_density_per_mb = round(total_variants / (genome_size_bp / 1_000_000), 3)
        snp_density_per_mb = round(snp_count / (genome_size_bp / 1_000_000), 3)
        indel_density_per_mb = round(indel_count / (genome_size_bp / 1_000_000), 3)
    else:
        variant_density_per_mb = "NA"
        snp_density_per_mb = "NA"
        indel_density_per_mb = "NA"

    # Per-chromosome Ti/Tv
    chr_titv = {}
    for chrom in chrom_counts.keys():
        chr_ti = chr_tv = 0
        with open(vcf_path, "r") as f2:
            for line in f2:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                c, _, _, ref, alt = parts[:5]
                if c != chrom:
                    continue
                alt_alleles = alt.split(",")
                if is_snp(ref, alt_alleles):
                    for a in alt_alleles:
                        if is_transition(ref, a):
                            chr_ti += 1
                        else:
                            chr_tv += 1
        chr_titv[chrom] = round(chr_ti / chr_tv, 3) if chr_tv > 0 else "NA"

    # Per-chromosome density
    chr_density = {}
    if chrom_sizes:
        for chrom, counts in chrom_counts.items():
            if chrom in chrom_sizes and chrom_sizes[chrom] > 0:
                size_mb = chrom_sizes[chrom] / 1_000_000
                chr_density[chrom] = {
                    "total_per_mb": round(counts["total"] / size_mb, 3),
                    "snps_per_mb": round(counts["snps"] / size_mb, 3),
                    "indels_per_mb": round(counts["indels"] / size_mb, 3),
                }
            else:
                chr_density[chrom] = {
                    "total_per_mb": "NA",
                    "snps_per_mb": "NA",
                    "indels_per_mb": "NA",
                }

    return {
        "total_variants": total_variants,
        "snp_count": snp_count,
        "indel_count": indel_count,
        "snp_prop": snp_prop,
        "indel_prop": indel_prop,
        "snp_indel_ratio": snp_indel_ratio,
        "titv_ratio": titv_ratio,
        "avg_qual": round(avg_qual, 3),
        "largest_indel": largest_indel,
        "mean_indel_size": round(mean_indel_size, 3),
        "insertions": insertion_count,
        "deletions": deletion_count,
        "multiallelic_count": multiallelic_count,
        "avg_snps_per_chr": round(avg_snps_per_chr, 2),
        "avg_indels_per_chr": round(avg_indels_per_chr, 2),
        "chrom_counts": dict(chrom_counts),
        "chr_titv": chr_titv,
        "mutation_spectrum": mutation_spectrum,
        "variant_density_per_mb": variant_density_per_mb,
        "snp_density_per_mb": snp_density_per_mb,
        "indel_density_per_mb": indel_density_per_mb,
        "chr_density": chr_density,
    }


def write_individual_report(sample_name, stats):
    path = os.path.join(OUTPUT_DIR, f"{sample_name}_summary.txt")
    with open(path, "w") as out:
        out.write(f"=== {sample_name} ===\n")
        out.write(f"Total variants: {stats['total_variants']}\n")
        out.write(f"Total SNPs: {stats['snp_count']} ({stats['snp_prop']}%)\n")
        out.write(f"Total INDELs: {stats['indel_count']} ({stats['indel_prop']}%)\n")
        out.write(f"SNP/INDEL ratio: {stats['snp_indel_ratio']}\n")
        out.write(f"Ti/Tv ratio: {stats['titv_ratio']}\n")
        out.write(f"Average QUAL: {stats['avg_qual']}\n")
        out.write(f"Largest INDEL size: {stats['largest_indel']}\n")
        out.write(f"Mean INDEL size: {stats['mean_indel_size']}\n")
        out.write(f"Insertions: {stats['insertions']}\n")
        out.write(f"Deletions: {stats['deletions']}\n")
        out.write(f"Multi-allelic sites: {stats['multiallelic_count']}\n")
        out.write(f"Avg SNPs per chromosome: {stats['avg_snps_per_chr']}\n")
        out.write(f"Avg INDELs per chromosome: {stats['avg_indels_per_chr']}\n")
        out.write(f"Variant density (per Mb): {stats['variant_density_per_mb']}\n")
        out.write(f"SNP density (per Mb): {stats['snp_density_per_mb']}\n")
        out.write(f"INDEL density (per Mb): {stats['indel_density_per_mb']}\n\n")

        out.write("Mutation Spectrum:\n")
        for mut_class, count in stats["mutation_spectrum"].items():
            out.write(f"  {mut_class}: {count}\n")

        out.write("\nPer Chromosome Counts:\n")
        for chrom, counts in stats["chrom_counts"].items():
            out.write(
                f"  {chrom}: Total={counts['total']}, SNPs={counts['snps']}, INDELs={counts['indels']}\n"
            )

        out.write("\nPer Chromosome Ti/Tv:\n")
        for chrom, val in stats["chr_titv"].items():
            out.write(f"  {chrom}: Ti/Tv={val}\n")

        if stats["chr_density"]:
            out.write("\nPer Chromosome Density (per Mb):\n")
            for chrom, vals in stats["chr_density"].items():
                out.write(
                    f"  {chrom}: Total={vals['total_per_mb']}, SNPs={vals['snps_per_mb']}, INDELs={vals['indels_per_mb']}\n"
                )

    print(f"Individual report saved: {path}")


def write_combined_summary(stats_dict):
    combined_path = os.path.join(OUTPUT_DIR, COMBINED_FILE)
    with open(combined_path, "w") as out:
        for sample_name, stats in stats_dict.items():
            out.write(f"=== {sample_name} ===\n")
            out.write(f"Total variants: {stats['total_variants']}\n")
            out.write(f"Total SNPs: {stats['snp_count']} ({stats['snp_prop']}%)\n")
            out.write(f"Total INDELs: {stats['indel_count']} ({stats['indel_prop']}%)\n")
            out.write(f"SNP/INDEL ratio: {stats['snp_indel_ratio']}\n")
            out.write(f"Ti/Tv ratio: {stats['titv_ratio']}\n")
            out.write(f"Average QUAL: {stats['avg_qual']}\n")
            out.write(f"Largest INDEL size: {stats['largest_indel']}\n")
            out.write(f"Mean INDEL size: {stats['mean_indel_size']}\n")
            out.write(f"Insertions: {stats['insertions']}\n")
            out.write(f"Deletions: {stats['deletions']}\n")
            out.write(f"Multi-allelic sites: {stats['multiallelic_count']}\n")
            out.write(f"Avg SNPs per chromosome: {stats['avg_snps_per_chr']}\n")
            out.write(f"Avg INDELs per chromosome: {stats['avg_indels_per_chr']}\n")
            out.write(f"Variant density (per Mb): {stats['variant_density_per_mb']}\n")
            out.write(f"SNP density (per Mb): {stats['snp_density_per_mb']}\n")
            out.write(f"INDEL density (per Mb): {stats['indel_density_per_mb']}\n\n")

            out.write("Mutation Spectrum:\n")
            for mut_class, count in stats["mutation_spectrum"].items():
                out.write(f"  {mut_class}: {count}\n")

            out.write("\nPer Chromosome Counts:\n")
            for chrom, counts in stats["chrom_counts"].items():
                out.write(
                    f"  {chrom}: Total={counts['total']}, SNPs={counts['snps']}, INDELs={counts['indels']}\n"
                )

            out.write("\nPer Chromosome Ti/Tv:\n")
            for chrom, val in stats["chr_titv"].items():
                out.write(f"  {chrom}: Ti/Tv={val}\n")

            if stats["chr_density"]:
                out.write("\nPer Chromosome Density (per Mb):\n")
                for chrom, vals in stats["chr_density"].items():
                    out.write(
                        f"  {chrom}: Total={vals['total_per_mb']}, SNPs={vals['snps_per_mb']}, INDELs={vals['indels_per_mb']}\n"
                    )

            out.write("\n\n")

    print(f"Combined summary saved: {combined_path}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    vcf_paths = load_vcf_paths()
    chrom_sizes = load_chrom_sizes(CHROM_SIZE_FILE)
    all_stats = {}

    if not vcf_paths:
        print("No VCF files found in INPUT_VCFS or INPUT_VCF_LIST_FILE.")
        return

    for vcf in vcf_paths:
        if not os.path.exists(vcf):
            print(f"File not found: {vcf}")
            continue

        sample_name = os.path.basename(vcf).replace(".vcf", "")
        print(f"Processing {sample_name}...")

        stats = parse_vcf(
            vcf_path=vcf,
            chrom_sizes=chrom_sizes if chrom_sizes else None,
            genome_size_bp=GENOME_SIZE_BP,
        )

        all_stats[sample_name] = stats
        write_individual_report(sample_name, stats)

    if all_stats:
        write_combined_summary(all_stats)
    else:
        print("No valid VCF files were processed.")


if __name__ == "__main__":
    main()
