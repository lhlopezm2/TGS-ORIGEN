#!/usr/bin/bash
#SBATCH -J 2-1-1
#SBATCH -D .
#SBATCH -o results/out-2-1-1.txt
#SBATCH -e results/err-2-1-1.txt
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --account=coffea_genomes
#SBATCH --time=1-23:00:00
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luish.lopezm@autonoma.edu.co
source ~/.bashrc

CPU=20
setup_v="-2-1-1"
fastq_v=${setup_v:0:2}
bam_v=${setup_v:0:4}
home="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis"
source "${home}/resource_usage.sh"
inter_med=5

fastq="base-calling/fastq${fastq_v}.fastq"
bam="alignment/bam${bam_v}.bam"
ref_fasta="GRCh38.fa"
vcf_folder="${home}/variant-calling/vcf${setup_v}"
vcf="${vcf_folder}/vcf${setup_v}.vcf.gz"
vcf_truth="HG002_GRCh38_benchmark.vcf.gz"
bed="GRCh38.bed"
metrics_prefix="metrics/metrics${setup_v}/metrics${setup_v}"

if [ -e "$vcf" ]; then
  echo "$vcf already exists."
else
  measure_resource_usage "${home}" "${inter_med}" "setup${setup_v}-variant_calling.txt" &
  measure_pid=$!
  echo "PID measure $measure_pid" >> "setup${setup_v}-variant_calling.txt" # Escribir el PID en el LOGFILE
  sleep "$((inter_med + 2))"

  echo "----------------"
  echo "Variant calling step"
  module load singularity
  module load nextflow/23.04.1
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Variant calling step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nextflow run wf-human-variation \
      -w ${vcf_folder} \
      -profile singularity \
      --snp --sv \
      --fast5_dir ./raw_reads\
      --bed ${bed} \
      --ref ${ref_fasta} \
      --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0'  \
      --sample_name "vcf${setup_v}" \
      --out_dir ${vcf_folder}
  module unload nextflow/23.04.1
  module unload singularity
  kill $measure_pid
fi

