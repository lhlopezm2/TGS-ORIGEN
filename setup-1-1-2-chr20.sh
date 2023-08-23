#!/usr/bin/bash
#SBATCH -J 1-1-2-chr20
#SBATCH -D .
#SBATCH -o results/out-1-1-2-chr20.txt
#SBATCH -e results/err-1-1-2-chr20.txt
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
vcf_folder="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/variant-calling/vcf-1-1-2-chr20"
vcf_prefix="vcf-1-1-2-chr20"
vcf="${vcf_folder}/${vcf_prefix}.vcf.gz"
vcf_truth="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/HG002_GRCh38_chr20.vcf.gz"
bed_truth="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/GRCh38_chr20.bed"
metrics_prefix="metrics/metrics-1-1-2-chr20/metrics-1-1-2-chr20"


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
HGREF=$ref_fasta /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compute metrics - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' hap.py $vcf_truth $vcf \
    --threads $CPU \
    -o $metrics_prefix \
    -T $bed_truth \
    -l "chr20"
conda deactivate

