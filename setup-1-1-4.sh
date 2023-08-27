#!/usr/bin/bash
#SBATCH -J 1-1-4
#SBATCH -D .
#SBATCH -o results/out-1-1-4.txt
#SBATCH -e results/err-1-1-4.txt
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
setup_v="-1-1-4"
fastq_v=${setup_v:0:2}
bam_v=${setup_v:0:4}
home="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis"
source "${home}/resource_usage.sh"
inter_med=5

fastq="base-calling/fastq${fastq_v}.fastq"
bam="alignment/bam${bam_v}.bam"
ref_fasta="GRCh38.fa"
vcf_folder="${home}/variant-calling/vcf${setup_v}"
vcf="${vcf_folder}/vcf${setup_v}.wf_snp.vcf.gz"
vcf_truth="HG002_GRCh38_benchmark.vcf.gz"
bed="GRCh38.bed"
metrics_prefix="metrics/metrics${setup_v}/metrics${setup_v}"

if [ -e "$fastq" ]; then
  echo "$fastq already exists."
else
  measure_resource_usage "${home}" "${inter_med}" "setup${setup_v}-base_calling.txt" &
  measure_pid=$!
  echo "PID measure $measure_pid" >> "setup${setup_v}-base_calling.txt" # Escribir el PID en el LOGFILE
  sleep "$((inter_med + 2))"

  echo "----------------"
  echo "Base-calling step"
  module load singularity
  module load guppy/6.4.6-gpu
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Base Calling - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' guppy_basecaller --disable_pings\
    -i ./raw_reads\
    -s "base-calling/fastq${fastq_v}"\
    --cpu_threads_per_caller 4\
    --flowcell FLO-PRO002M\
    --kit SQK-RBK112-96\
    --recursive -x 'cuda:all:50G'\
    --num_callers 5\
    --compress_fastq 
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Merge fastq files - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' zcat "base-calling/fastq${fastq_v}/pass/fastq_runid*.fastq.gz" > $fastq
  module unload guppy/6.4.6-gpu
  module unload singularity
  kill $measure_pid
fi

if [ -e "$bam" ]; then
  echo "$bam already exists."
else
  measure_resource_usage "${home}" "${inter_med}" "setup${setup_v}-alignment.txt" &
  measure_pid=$!
  echo "PID measure $measure_pid" >> "setup${setup_v}-alignment.txt" # Escribir el PID en el LOGFILE
  sleep "$((inter_med + 2))"

  echo "----------------"
  echo "Alignment step"
  module load minimap2/2.24
  module load samtools/1.15.1
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Alignment step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' minimap2 -a -z 600,200 -x map-ont $ref_fasta $fastq -t $CPU \
    |samtools view -Shu |samtools sort -@ $CPU -o $bam --output-fmt BAM
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Bam indexing - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' samtools index $bam -@ $CPU
  module unload minimap2/2.24
  module unload samtools/1.15.1
  kill $measure_pid
fi

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
      --bam ${bam} \
      --bed ${bed} \
      --ref ${ref_fasta} \
      --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0'  \
      --sample_name "vcf${setup_v}" \
      --out_dir ${vcf_folder}
  module unload nextflow/23.04.1
  module unload singularity
  kill $measure_pid
fi

echo "----------------"
echo "Compute metrics"
conda activate happy
HGREF=$ref_fasta /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compute metrics - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' hap.py $vcf_truth $vcf \
    --threads $CPU \
    -o $metrics_prefix \
    -T $bed
conda deactivate

