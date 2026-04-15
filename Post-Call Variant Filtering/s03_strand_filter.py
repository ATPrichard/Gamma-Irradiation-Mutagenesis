#!/usr/bin/env python3

import os
import sys

# CONFIG SECTION

INPUT_LIST = "VCF_lists/HF_vcf_list.txt"
RESULTS_DIR = "Results"

MIN_ADF = 2
MIN_ADR = 2
PROGRESS_UPDATE_INTERVAL = 50



def parse_sample_field(sample_field):
    if sample_field is None or sample_field == ".":
        return []
    return sample_field.split(":")


def parse_two_allele_counts(raw):
    if raw is None or raw == "." or raw == "":
        return (-1, -1)
    parts = raw.split(",")
    if len(parts) < 2:
        return (-1, -1)
    try:
        return int(parts[0]), int(parts[1])
    except Exception:
        return (-1, -1)


def safe_basename_noext(path):
    base = os.path.basename(path)
    if base.lower().endswith(".vcf"):
        return base[:-4]
    return os.path.splitext(base)[0]


def process_single_vcf(vcf_path, soi_line_index):
    if not os.path.isfile(vcf_path):
        print(f"WARNING: input VCF not found, skipping: {vcf_path}", file=sys.stderr)
        return None  # Return None to indicate skipping

    with open(vcf_path, "r") as fh:
        headers = []
        samples = []

        for line in fh:
            if line.startswith("##"):
                headers.append(line)
            elif line.startswith("#CHROM"):
                headers.append(line)
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                break

        variant_lines = fh.readlines()

    num_samples = len(samples)
    soi_col_index = 9 + soi_line_index
    soi_sample_idx = soi_col_index - 10

    if soi_sample_idx < 0 or soi_sample_idx >= num_samples:
        print(f"ERROR: Invalid SOI index for {vcf_path}", file=sys.stderr)
        return None

    sample_name = samples[soi_sample_idx]
    base_noext = safe_basename_noext(vcf_path)

    out_vcf = os.path.join(RESULTS_DIR, f"{base_noext}_SF.vcf")
    summary_file = os.path.join(RESULTS_DIR, f"{base_noext}_SF_summary.txt")

    print(f"Processing: {vcf_path}")
    print(f"  SOI: {sample_name}")
    print(f"  Output: {out_vcf}")

    total = kept = removed = 0

    with open(out_vcf, "w") as out:
        for h in headers:
            out.write(h)

        for i, line in enumerate(variant_lines, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                removed += 1
                total += 1
                continue

            format_field = parts[8]
            sample_fields = parts[9:]

            format_keys = format_field.split(":")
            fmt_index = {k: idx for idx, k in enumerate(format_keys)}

            if soi_sample_idx >= len(sample_fields):
                removed += 1
                total += 1
                continue

            soi_values = parse_sample_field(sample_fields[soi_sample_idx])

            # Extract ADF / ADR
            adf_ref = adf_alt = -1
            adr_ref = adr_alt = -1

            if "ADF" in fmt_index:
                idx = fmt_index["ADF"]
                if idx < len(soi_values):
                    adf_ref, adf_alt = parse_two_allele_counts(soi_values[idx])

            if "ADR" in fmt_index:
                idx = fmt_index["ADR"]
                if idx < len(soi_values):
                    adr_ref, adr_alt = parse_two_allele_counts(soi_values[idx])

            # Apply strand filter
            pass_filters = (
                adf_ref >= MIN_ADF and
                adf_alt >= MIN_ADF and
                adr_ref >= MIN_ADR and
                adr_alt >= MIN_ADR
            )

            total += 1

            if pass_filters:
                out.write(line)
                kept += 1
            else:
                removed += 1

            if i % PROGRESS_UPDATE_INTERVAL == 0:
                print(f"  Processed {i:,} variants...")

    with open(summary_file, "w") as s:
        s.write(f"Input VCF: {vcf_path}\n")
        s.write(f"Output VCF: {out_vcf}\n")
        s.write(f"Sample of Interest (SOI): {sample_name}\n\n")

        s.write("Filtering criteria:\n")
        s.write("  Variant type: SNPs only (pipeline-split upstream)\n")
        s.write("  Strand balance requirements (SOI):\n")
        s.write(f"    - ADF (REF) >= {MIN_ADF}\n")
        s.write(f"    - ADF (ALT) >= {MIN_ADF}\n")
        s.write(f"    - ADR (REF) >= {MIN_ADR}\n")
        s.write(f"    - ADR (ALT) >= {MIN_ADR}\n")
        s.write("  Missing or malformed ADF/ADR values: FAIL\n\n")

        s.write(f"Total sites processed: {total}\n")
        s.write(f"Sites kept: {kept}\n")
        s.write(f"Sites discarded: {removed}\n")

        if total > 0:
            s.write(f"Percent kept: {kept / total * 100:.2f}%\n")
            s.write(f"Percent discarded: {removed / total * 100:.2f}%\n")

    print(f"Finished {sample_name}: kept {kept}, discarded {removed}\n")

    return {"total": total, "kept": kept, "removed": removed}


def main():
    if not os.path.isfile(INPUT_LIST):
        print(f"ERROR: INPUT_LIST not found: {INPUT_LIST}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(RESULTS_DIR, exist_ok=True)

    global_stats = {"total": 0, "kept": 0, "removed": 0}
    sample_summaries = []

    with open(INPUT_LIST, "r") as f:
        lines = [ln.rstrip("\n") for ln in f.readlines()]

    line_num = 0
    for raw in lines:
        if not raw or raw.startswith("#"):
            continue
        line_num += 1
        stats = process_single_vcf(raw.strip(), line_num)
        if stats:
            sample_summaries.append(stats)
            for key in global_stats:
                global_stats[key] += stats[key]

    global_summary_file = os.path.join(RESULTS_DIR, "Global_SF_summary.txt")
    with open(global_summary_file, "w") as g:
        g.write("Global Strand Filter Summary\n")
        g.write("============================\n\n")
        g.write(f"Total sites processed across all samples: {global_stats['total']}\n")
        g.write(f"Total sites kept across all samples: {global_stats['kept']}\n")
        g.write(f"Total sites discarded across all samples: {global_stats['removed']}\n")

        if global_stats['total'] > 0:
            percent_kept_list = [(s["kept"] / s["total"] * 100) for s in sample_summaries]
            mean_percent_kept = sum(percent_kept_list) / len(percent_kept_list)
            g.write(f"Mean percent kept across samples: {mean_percent_kept:.2f}%\n")

    print("All files processed. Global summary written.")


if __name__ == "__main__":
    main()
