#!/usr/bin/env python3

import os

# CONFIG SECTION

INPUT_VCF = "Results/raw_GT.vcf"

RESULTS_DIR = "Results"
SUMMARY_DIR = "Summaries"

WILDTYPE_IDX = # 0 based, tell where wildtype sample is, if none then comment out

PROGRESS_UPDATE_INTERVAL = 100000

HOMOZYGOUS_THRESHOLD = 0.9


def parse_genotype(field):
    return field.split(":")[0] if ":" in field else field


def is_heterozygous(gt):
    return gt in ("0/1", "1/0")



def process_vcf(input_vcf, results_dir, summary_dir, wt_idx,
                progress_interval, homo_threshold):

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(summary_dir, exist_ok=True)

    global_summary_data = []

    with open(input_vcf, "r") as infile:
        headers = []
        samples = []

        for line in infile:
            if line.startswith("##"):
                headers.append(line)
            elif line.startswith("#CHROM"):
                headers.append(line)
                header_parts = line.strip().split("\t")
                samples = header_parts[9:]
                break

        treated_samples = [i for i in range(len(samples)) if i != wt_idx]
        variant_lines = infile.readlines()

    print(f"Detected {len(samples)} samples total.")
    print(f"Homozygous agreement threshold: {homo_threshold * 100:.1f}%")
    print("-" * 60)


    for t_idx in treated_samples:

        treated_name = samples[t_idx]
        sample_clean = os.path.basename(treated_name).replace(".", "_")

        output_vcf = os.path.join(results_dir, f"{sample_clean}_HF.vcf")
        summary_file = os.path.join(summary_dir, f"{sample_clean}_HF_summary.txt")

        kept = 0
        discarded = 0
        total = 0

        print(f"Processing treated sample: {treated_name}")

        with open(output_vcf, "w") as out:
            for h in headers:
                out.write(h)

            for i, line in enumerate(variant_lines, start=1):

                if not line.strip() or line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                sample_fields = parts[9:]
                gts = [parse_genotype(f) for f in sample_fields]

                treated_gt = gts[t_idx]
                total += 1

                # SOI must be heterozygous
                if not is_heterozygous(treated_gt):
                    discarded += 1
                    continue

                # Collect non-SOI genotypes
                non_soi_gts = [gt for j, gt in enumerate(gts) if j != t_idx]
                non_soi_gts = [gt for gt in non_soi_gts if gt != "./."]

                if len(non_soi_gts) == 0:
                    discarded += 1
                    continue

                count_00 = non_soi_gts.count("0/0")
                count_11 = non_soi_gts.count("1/1")

                max_hom = max(count_00, count_11)
                proportion = max_hom / len(non_soi_gts)

                if proportion >= homo_threshold:
                    out.write(line)
                    kept += 1
                else:
                    discarded += 1

                if i % progress_interval == 0:
                    print(f"  Processed {i:,} lines...")

        percent_kept = (kept / total * 100) if total > 0 else 0

        with open(summary_file, "w") as s:
            s.write(f"Input VCF: {input_vcf}\n")
            s.write(f"Homozygous agreement threshold: {homo_threshold * 100:.1f}%\n\n")
            s.write(f"Total sites processed: {total}\n")
            s.write(f"Sites kept: {kept}\n")
            s.write(f"Sites discarded: {discarded}\n")
            s.write(f"Percent kept: {percent_kept:.2f}%\n")

        global_summary_data.append({
            "sample": treated_name,
            "total": total,
            "kept": kept,
            "discarded": discarded,
            "percent_kept": percent_kept
        })

        print(f"Finished {treated_name}: kept {kept}, discarded {discarded}")
        print("-" * 60)


    global_summary_file = os.path.join(summary_dir, "HF_global_summary.txt")

    percent_values = [entry["percent_kept"] for entry in global_summary_data]
    kept_values = [entry["kept"] for entry in global_summary_data]

    if percent_values:
        mean_percent = sum(percent_values) / len(percent_values)

        sorted_percents = sorted(percent_values)
        mid = len(sorted_percents) // 2
        if len(sorted_percents) % 2 == 0:
            median_percent = (sorted_percents[mid - 1] + sorted_percents[mid]) / 2
        else:
            median_percent = sorted_percents[mid]

        min_percent = min(percent_values)
        max_percent = max(percent_values)
        mean_kept = sum(kept_values) / len(kept_values)
    else:
        mean_percent = median_percent = min_percent = max_percent = mean_kept = 0

    with open(global_summary_file, "w") as g:
        g.write("VCF Heterozygous Site Finder - Global Summary\n")
        g.write("=" * 60 + "\n\n")
        g.write(f"Input VCF: {input_vcf}\n")
        g.write(f"Homozygous agreement threshold: {homo_threshold * 100:.1f}%\n")
        g.write(f"Total samples processed: {len(global_summary_data)}\n\n")

        g.write("Per-Sample Results\n")
        g.write("-" * 60 + "\n")

        for entry in global_summary_data:
            g.write(f"Sample: {entry['sample']}\n")
            g.write(f"  Total sites: {entry['total']}\n")
            g.write(f"  Sites kept: {entry['kept']}\n")
            g.write(f"  Sites discarded: {entry['discarded']}\n")
            g.write(f"  Percent kept: {entry['percent_kept']:.2f}%\n")
            g.write("-" * 40 + "\n")

        g.write("\nOverall Statistics\n")
        g.write("=" * 60 + "\n")
        g.write(f"Mean percent kept: {mean_percent:.2f}%\n")
        g.write(f"Median percent kept: {median_percent:.2f}%\n")
        g.write(f"Minimum percent kept: {min_percent:.2f}%\n")
        g.write(f"Maximum percent kept: {max_percent:.2f}%\n")
        g.write(f"Mean sites kept (absolute): {mean_kept:.2f}\n")

    print("All treated samples processed successfully!")
    print(f"Global summary written to: {global_summary_file}")



if __name__ == "__main__":
    process_vcf(
        INPUT_VCF,
        RESULTS_DIR,
        SUMMARY_DIR,
        WILDTYPE_IDX,
        PROGRESS_UPDATE_INTERVAL,
        HOMOZYGOUS_THRESHOLD
    )
