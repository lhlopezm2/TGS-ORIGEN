#!/usr/bin/bash
#SBATCH -J 1-1-4-chr20
#SBATCH -D .
#SBATCH -o results/out-1-1-4-chr20.txt
#SBATCH -e results/err-1-1-4-chr20.txt
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
bam="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/HG002_GRCh38_chr20.bam"
ref_fasta="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/GRCh38_chr20.fa"
vcf_folder="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/variant-calling/vcf-1-1-4-chr20"
vcf="${vcf_folder}/vcf-1-1-4-chr20.wf_snp.vcf.gz"
vcf_truth="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/HG002_GRCh38_chr20.vcf.gz"
bed_truth="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/GRCh38_chr20.bed"
metrics_prefix="metrics/metrics-1-1-4-chr20/metrics-1-1-4-chr20"

if [ -e "$vcf" ]; then
  echo "$vcf already exists."
else
  echo "----------------"
  echo "Variant calling step"
  module load singularity
  module load nextflow/23.04.1
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Variant calling step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nextflow run wf-human-variation \
      -w ${vcf_folder} \
      -profile singularity \
      --snp --sv \
      --bam ${bam} \
      --bed ${bed_truth} \
      --ref ${ref_fasta} \
      --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0'  \
      --sample_name vcf-1-1-4-chr20 \
      --out_dir ${vcf_folder}
  module unload nextflow/23.04.1
  module unload singularity
fi

echo "----------------"
echo "Compute metrics"
conda activate happy
HGREF=$ref_fasta /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compute metrics - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' hap.py $vcf_truth $vcf \
    --threads $CPU \
    -o $metrics_prefix \
    -T $bed_truth \
    -l "chr20"
conda deactivate
