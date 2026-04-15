# Bioinformatic Pipeline for the Analysis of Gamma Irradiated *Daphnia*
---
**Author**: Andrew Prichard \
**Affiliation**: University of Missouri - Columbia, Department of Biological Sciences \

#### As part of the thesis project entitled "Characterization of the rate and spectrum of gamma radiation mutagenesis in the microcrustacean *Daphnia*"
---

## Variant Calling Pipeline:
Follows a standard genomic workflow, code is included in the "Variant Calling Pipeline" directory.

**Follows these basic steps:** 
1. Index reference genome (only needs to be done once with each genome)
1. Quality control of reads (Fastp Trim)
2. Align reads to reference genome (BWA-Mem, Bowtie2)
3. Read processing (Samtools; fixmate, sort, markdup)
4. Variant calling (BCFTools; mpileup, call)
5. Preliminary filtering (BCFTools; filter)

---

## Post-Call Variant Filtering:
- Custom solution for analyzing the *Daphnia* genome after mutagenesis
- Script names are included in their relevant info section below. Scripts are named with a prefix in their folder as; s0X_ = SNP script, i0X_ = INDEL script, b0X_ = Both

### Step 1.) Pre-Processing:
Accepts a raw VCF file containing variant data generated at the end of the "Variant Calling Pipeline". What this script does is separates the original VCF file into two separate files, one containing SNP calls and the other containing INDEL calls. 

#### Scripts associated with this step:
- split.py 
---
### Step 2.) Purity Scripts:
Ensures proper, high confidence genotype calls. Looks at homozygous genotype (GT) calls and ensures complete homozygosity (i.e. no reads supporting alternate allele, if called for REF all reads should support REF and not ALT). If a homozygous GT that has alternate allele support is identified, then that GT will be masked. If X% of samples have a masked genotype at a given position, then the whole site is removed. (For our analysis we used a threshold of 70%).

#### Scripts associated with this step:
- genotype_check.py (SNP)
- homozygous_consistency.py (INDEL)
---

 ### Step 3.) Comparison Scripts:
 Iterates through all samples present in the file (-WT) and defines the sample as the sample of interest (SOI). Then looks for sites where the SOI is called heterozygous. The SOI is then compared to all other samples present in the file and X% (we used 90%) of sites need to be homozygous, if the site does not fulfill that requirement, then the site is thrown out and moves on to the next. Generates a VCF for each sample.

 #### Scripts associated with this step:
 - heterozygous_filter.py (SNP)
 - comparison_script.py (INDEL)
---
### Step 4.) Merge:
Takes the SNP file for a given sample and INDEL file for a sample and combines them together in position order.

 #### Scripts associated with this step:
 - combine.py
---
### Step 5.) Basic Statistic Generation:
Generates useful stats for each of the files (i.e. number of total variants, SNP and INDEL count, Ti/Tv, per chromosome breakdowns, plus a lot more).
 #### Scripts associated with this step:
 - stat_summary.py
---

