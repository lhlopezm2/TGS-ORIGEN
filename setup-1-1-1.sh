#!/usr/bin/bash
#SBATCH -J 1-1-1
#SBATCH -D .
#SBATCH -o results/out-1-1-1.txt
#SBATCH -e results/err-1-1-1.txt
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
setup_v="-1-1-1"
fastq_v=${setup_v:0:2}
bam_v=${setup_v:0:4}
home="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis"


fastq="base-calling/fastq${fastq_v}.fastq"
bam="alignment/bam${bam_v}.bam"
ref_fasta="GRCh38.fa"
MODEL_NAME="r941_prom_hac_g360+g422"
vcf_folder="${home}/variant-calling/vcf${setup_v}"
vcf="${vcf_folder}/merge_output.vcf.gz"
vcf_truth="HG002_GRCh38_benchmark.vcf.gz"
bed="GRCh38.bed"
metrics_prefix="metrics/metrics${setup_v}/metrics${setup_v}"

if [ -e "$fastq" ]; then
  echo "$fastq already exists."
else
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
  
fi

if [ -e "$bam" ]; then
  echo "$bam already exists."
else
  echo "----------------"
  echo "Alignment step"
  module load minimap2/2.24
  module load samtools/1.15.1
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Alignment step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' minimap2 -a -z 600,200 -x map-ont $ref_fasta $fastq -t $CPU \
    |samtools view -Shu |samtools sort -@ $CPU -o $bam --output-fmt BAM
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Bam indexing - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' samtools index $bam -@ $CPU
  module unload minimap2/2.24
  module unload samtools/1.15.1
fi

if [ -e "$vcf" ]; then
  echo "$vcf already exists."
else
  echo "----------------"
  echo "Variant calling step"
  conda activate clair3
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Variant calling step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' run_clair3.sh \
    --bam_fn=$bam \
    --ref_fn=$ref_fasta \
    --threads=$CPU \
    --platform="ont" \
    --model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \
    --output=$vcf_folder
  conda deactivate
fi

echo "----------------"
echo "Compute metrics"
conda activate happy
HGREF=$ref_fasta /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compute metrics - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' hap.py $vcf_truth $vcf \
    --threads $CPU \
    -o $metrics_prefix \
    -T $bed
conda deactivate

