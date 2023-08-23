#!/usr/bin/bash
#SBATCH -J 1-1-1-chr20
#SBATCH -D .
#SBATCH -o out-1-1-1-chr20.txt
#SBATCH -e err-1-1-1-chr20.txt
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
MODEL_NAME="r941_prom_hac_g360+g422"
vcf_folder="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/variant-calling/vcf-1-1-1-chr20"
vcf="${vcf_folder}/merge_output.vcf.gz"
vcf_truth="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/HG002_GRCh38_chr20.vcf.gz"
bed_truth="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/HG002_GRCh38_chr20.bed"
metrics_prefix="metrics/metrics-1-1-1-chr20"

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
    -T $bed_truth \
    -l "chr20"
conda deactivate

