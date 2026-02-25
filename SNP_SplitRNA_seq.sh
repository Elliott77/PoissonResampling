#!/bin/bash
## SNP_SplitRNA_seq.sh
## Elliott Ferris, 2019
##
## Split allele-specific RNA-seq BAM files using SNPsplit, then sort and index.
## Requires: SNPsplit, samtools
##
## Usage:
##   1. Edit the CONFIGURATION section below.
##   2. Submit via SLURM:  sbatch SNP_SplitRNA_seq.sh
##      Or run directly:   bash SNP_SplitRNA_seq.sh

###############################################################################
## SLURM OPTIONS — modify for your cluster
###############################################################################
#SBATCH --time=6:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=<your-account>
#SBATCH --job-name=snpsplit_rna
#SBATCH -o logs/split_rna-%j.out-%N
#SBATCH -e logs/split_stderr-%j.txt

###############################################################################
## CONFIGURATION — edit these variables for your experiment
###############################################################################

## Directory containing input BAM files
BAM_DIR="<path to BAM files>"

## Path to SNPsplit executable
SNPSPLIT="<path to SNPsplit>/SNPsplit"

## Path to strain-specific SNP file (e.g., from the Mouse Genomes Project)
## See: https://www.sanger.ac.uk/data/mouse-genomes-project/
SNP_FILE="<path to SNP file>/all_SNPs_CAST_EiJ_GRCm38.txt.gz"

## Number of threads for samtools sort
THREADS=24

###############################################################################
## SETUP
###############################################################################

set -euo pipefail

module load samtools 2>/dev/null || true
SAMTOOLS_PATH=$(which samtools)

cd "${BAM_DIR}"

###############################################################################
## SPLIT, SORT, AND INDEX EACH BAM FILE
###############################################################################

for bam in *.bam; do
    sample="${bam%.bam}"
    echo "Processing: ${sample}"

    ## Skip if already processed
    if [ -f "${sample}.allele_flagged.bam" ]; then
        echo "  Already split — skipping."
        continue
    fi

    ## Split reads by allele with SNPsplit
    echo "  Starting SNPsplit at $(date)"
    "${SNPSPLIT}" \
        --samtools_path "${SAMTOOLS_PATH}" \
        --conflicting \
        --paired \
        --snp_file "${SNP_FILE}" \
        "${bam}"

    ## Sort allele-specific BAMs
    echo "  Sorting genome1..."
    samtools sort -@ "${THREADS}" -o "${sample}_sort.genome1.bam" "${sample}.genome1.bam"
    echo "  Sorting genome2..."
    samtools sort -@ "${THREADS}" -o "${sample}_sort.genome2.bam" "${sample}.genome2.bam"

    ## Index sorted BAMs
    samtools index "${sample}_sort.genome1.bam"
    samtools index "${sample}_sort.genome2.bam"

    echo "  Done: ${sample}"
done

echo "All samples complete at $(date)"
