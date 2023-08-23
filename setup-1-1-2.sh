#!/usr/bin/bash
#SBATCH -J 1-1-2
#SBATCH -D .
#SBATCH -o results/out-1-1-2.txt
#SBATCH -e results/err-1-1-2.txt
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
fastq="base-calling/fastq-1.fastq"
bam="alignment/bam-1-1.bam"
ref_fasta="GRCh38.fa"
MODEL_NAME="r941_prom_hac_g360+g422"
vcf_folder="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/variant-calling/vcf-1-1-2"
vcf_prefix="vcf-1-1-2"
vcf="${vcf_folder}/${vcf_prefix}.vcf.gz"
vcf_truth="HG002_GRCh38_benchmark.vcf.gz"
bed_truth="GRCh38.bed"
bed_intersection="alignment/bed-intersection.bed"
bed_alignment="alignment/bed-alignment.bed"
metrics_prefix="metrics/metrics-1-1-2/metrics-1-1-2"
focus=0

if [ -e "$fastq" ]; then
  echo "$fastq already exists."
else
  echo "----------------"
  echo "Base-calling step"
  module load singularity
  module load guppy/6.4.6-gpu
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Base Calling - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' guppy_basecaller --disable_pings\
    -i ./raw_reads\
    -s base-calling/fastq-1\
    --cpu_threads_per_caller 4\
    --flowcell FLO-PRO002M\
    --kit SQK-RBK112-96\
    --recursive -x 'cuda:all:50G'\
    --num_callers 5\
    --compress_fastq 
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Merge fastq files - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' zcat base-calling/fastq-1/pass/fastq_runid*.fastq.gz > $fastq
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

if [ $focus -eq 1 ]; then
  if [ -e "$bed_intersection" ]; then
    echo "$bed_intersection already exists."
  else
    echo "----------------"
    echo "Intersect bam coverage with groundtruth bed file"
    conda activate clair3
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Bam to Bed - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bedtools bamtobed -i $bam > $bed_alignment
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Bed Intersection - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bedtools intersect -a $bed_alignment -b $bed_truth > $bed_intersection
    conda deactivate
  fi
else
  bed_intersection=$bed_truth
fi


if [ -e "$vcf" ]; then
  echo "$vcf already exists."
else
  echo "----------------"
  echo "Variant calling step"
  module load singularity
  module load tensorflow-gpu/2.6.2
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Variant calling step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' singularity exec --bind /usr/lib/locale/ \
    pepper_deepvariant_r0.8.sif \
    run_pepper_margin_deepvariant call_variant \
    -b "${bam}" \
    -f "${ref_fasta}" \
    -o "${vcf_folder}" \
    -p "${vcf_prefix}" \
    -t "${CPU}" \
    --ont_r9_guppy5_sup
  module unload tensorflow-gpu/2.6.2
  module unload singularity
fi

echo "----------------"
echo "Compute metrics"
conda activate happy
HGREF=$ref_fasta /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compute metrics - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' hap.py $vcf_truth $vcf --threads $CPU -o $metrics_prefix -T $bed_intersection
conda deactivate

