#!/usr/bin/bash
#SBATCH -J 1-1-3
#SBATCH -D .
#SBATCH -o results/out-1-1-3.txt
#SBATCH -e results/err-1-1-3.txt
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
fastq_index="base-calling/fastq-1.fastq.index"
bam="alignment/bam-1-1.bam"
ref_fasta="GRCh38.fa"
MODEL_NAME="r941_prom_hac_g360+g422"
vcf="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/variant-calling/vcf-1-1-3/vcf-1-1-3.vcf"
chr="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis/variant-calling/vcf-1-1-3"
#vcf="${vcf_folder}/merge_output.vcf.gz"
vcf_truth="HG002_GRCh38_benchmark.vcf.gz"
bed_truth="GRCh38.bed"
bed_intersection="alignment/bed-intersection.bed"
bed_alignment="alignment/bed-alignment.bed"
metrics_prefix="metrics/metrics-1-1-3/metrics-1-1-3"
focus=0
home="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis"

contigs=("chr1:1-248956422" "chr2:1-242193529" "chr3:1-198295559" "chr4:1-190214555" "chr5:1-181538259" "chr6:1-170805979" "chr7:1-159345973" "chr8:1-145138636" "chr9:1-138394717" "chr10:1-133797422" "chr11:1-135086622" "chr12:1-133275309" "chr13:1-114364328" "chr14:1-107043718" "chr15:1-101991189" "chr16:1-90338345" "chr17:1-83257441" "chr18:1-80373285" "chr19:1-58617616" "chr20:1-64444167" "chr21:1-46709983" "chr22:1-50818468")

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

if [ -e "$fastq_index" ]; then
  echo "$fastq_index already exists."
else
  echo "----------------"
  echo "Fastq indexing step"
  conda activate nanopolish
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Fastq indexing step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nanopolish index -d ./raw_reads $fastq
  conda deactivate
fi

if [ -e "${chr}/chr1.vcf" ]; then
  echo "chr1.vcf already exists."
else
  echo "----------------"
  echo "Variant calling step"
  conda activate nanopolish
  for region in "${contigs[@]}"; do
    contig=$(echo "$region" | cut -d':' -f1)
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Variant calling step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nanopolish variants \
      --outfile="${chr}/${contig}.vcf" --threads=$CPU\
      --ploidy=2\
      --reads=$fastq\
      --bam=$bam\
      --genome=$ref_fasta\
      --window="${region}"
  done
  conda deactivate
fi

if [ -e "${chr}/chr1.vcf.gz.tbi" ]; then
  echo "chr1.vcf.gz.tbi already exists."
else
  echo "----------------"
  echo "Compressing and indexing vcf files step"
  conda activate nanopolish
  for region in "${contigs[@]}"; do
    contig=$(echo "$region" | cut -d':' -f1)
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compressing vcf files step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bgzip -k -f "${chr}/${contig}.vcf"
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Indexing vcf files step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' tabix -p vcf "${chr}/${contig}.vcf.gz"
  done
  conda deactivate
fi

if [ -e "${vcf}" ]; then
  echo "${vcf} already exists."
else
  echo "----------------"
  echo "Merging vcf files"
  conda activate nanopolish
  cd $chr
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Merging vcf files step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bcftools merge chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz --force-samples -o $vcf
  cd $home
  conda deactivate
fi

conda activate nanopolish
/shared/home/sorozcoarias/anaconda3/bin/time -f 'Compressing merged vcf step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bgzip -k -f "${vcf}"
/shared/home/sorozcoarias/anaconda3/bin/time -f 'Indexing merged vcf step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' tabix -p vcf "${vcf}.gz"
conda deactivate

echo "----------------"
echo "Compute metrics"
conda activate happy
HGREF=$ref_fasta /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compute metrics - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' hap.py $vcf_truth $vcf --threads $CPU -o $metrics_prefix -T $bed_intersection
conda deactivate

