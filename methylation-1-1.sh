#!/usr/bin/bash
#SBATCH -J meth-1-1
#SBATCH -D .
#SBATCH -o results/out-meth-1-1.txt
#SBATCH -e results/err-meth-1-1.txt
#SBATCH -n 22
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
setup_v="-1-1-1"
fastq_v=${setup_v:0:2}
bam_v=${setup_v:0:4}
home="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis"
source "${home}/resource_usage.sh"
inter_med=5

fastq="base-calling/fastq${fastq_v}.fastq"
bam="alignment/bam${bam_v}.bam"
ref_fasta="GRCh38.fa"
tsv="methylation-calling/tsv${bam_v}.tsv"
fastq_index="base-calling/fastq${fastq_v}.fastq.index"

measure_resource_usage "${home}" "${inter_med}" "methylation${bam_v}.txt" &
measure_pid=$!
echo "PID measure $measure_pid" >> "methylation${bam_v}.txt" # Escribir el PID en el LOGFILE
sleep "$((inter_med + 2))"

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

kill $measure_pid
