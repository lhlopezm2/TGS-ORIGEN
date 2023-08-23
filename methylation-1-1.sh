#!/usr/bin/bash
#SBATCH -J 1-1-3
#SBATCH -D .
#SBATCH -o results/out-meth-1-1.txt
#SBATCH -e results/err-meth-1-1.txt
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --account=coffea_genomes
#SBATCH --time=1-23:00:00
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luish.lopezm@autonoma.edu.co
source ~/.bashrc

CPU=20
fastq="base-calling/fastq-1.fastq"
bam="alignment/bam-1-1.bam"
ref_fasta="GRCh38.fa"
tsv="methylation-calling/tsv-1-1.tsv"
fastq_index="base-calling/fastq-1.fastq.index"

if [ -e "$fastq_index" ]; then
  echo "$fastq_index already exists."
else
  echo "----------------"
  echo "Fastq indexing step"
  conda activate nanopolish
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Fastq indexing step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nanopolish index -d ./raw_reads $fastq
  conda deactivate
fi


if [ -e "${tsv}" ]; then
  echo "${tsv} already exists."
else
  echo "----------------"
  echo "Calling methylations"
  conda activate nanopolish
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Variant calling step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nanopolish call-methylation \
    --threads=$CPU\
    --reads=$fastq\
    --bam=$bam\
    --genome=$ref_fasta\
    --progress > ${tsv}
  conda deactivate
fi



