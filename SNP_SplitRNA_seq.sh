#!/bin/sh
## SNP_SplitRNA_seq.sh
## Elliott Ferris
## 01/05/2019
#SBATCH --time=6:30:00 # Walltime
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1       # We only run one R instance = 1 task
#SBATCH --account=<account>
#SBATCH --account=gregg-kp

#SBATCH --job-name=<Job Name>
#SBATCH -o <path to log directory>/split_rna-%j.out-%N
#SBATCH -e  <path to log directory>/split_stderr-%j.txt

module load samtools
STP=$(which samtools)
scratch_sam=<path to BAM files>
threads=24
counts=<path to counts>/tissue1_allele_counts.txt

cd $scratch_sam
Count=1
cd $scratch_sam
for f in `ls -1 *.bam | sed 's/.bam//'`
do
echo ${f}

if [ ! -f "${f}.allele_flagged.bam" ]; then
	## Sort Reads
	echo "Starting SNPsplit at `date`"
	cd $scratch_sam
	/uufs/chpc.utah.edu/common/home/u0368716/bin/SNPsplit_v0.3.2/SNPsplit --samtools_path $STP --conflicting --paired --snp_file <path to SNP file directory>/all_SNPs_CAST_EiJ_GRCm38.txt.gz ${f}.bam
	cd $scratch_sam
	samtools sort -@ $threads -o ${f}_sort.genome1.bam ${f}.genome1.bam
	samtools sort -@ $threads -o ${f}_sort.genome2.bam ${f}.genome2.bam
	samtools index ${f}_sort.genome1.bam
	samtools index ${f}_sort.genome2.bam
fi

echo done
done

echo "End of program at `date`"


