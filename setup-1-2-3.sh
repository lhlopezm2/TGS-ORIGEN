#!/usr/bin/bash
#SBATCH -J 1-2-3
#SBATCH -D .
#SBATCH -o results/out-1-2-3.txt
#SBATCH -e results/err-1-2-3.txt
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
setup_v="-1-2-3"
fastq_v=${setup_v:0:2}
bam_v=${setup_v:0:4}
home="/shared/home/sorozcoarias/coffea_genomes/Simon/Luis"
source "${home}/resource_usage.sh"
inter_med=5

fastq="base-calling/fastq${fastq_v}.fastq"
fastq_index="base-calling/fastq${fastq_v}.fastq.index"
bam="alignment/bam${bam_v}.bam"
ref_fasta="GRCh38.fa"
vcf="${home}/variant-calling/vcf${setup_v}/vcf${setup_v}.vcf"
chr="${home}/variant-calling/vcf${setup_v}"
vcf_truth="HG002_GRCh38_benchmark.vcf.gz"
bed="GRCh38.bed"
metrics_prefix="metrics/metrics${setup_v}/metrics${setup_v}"

contigs=("chr1:1-248956422" "chr2:1-242193529" "chr3:1-198295559" "chr4:1-190214555" "chr5:1-181538259" "chr6:1-170805979" "chr7:1-159345973" "chr8:1-145138636" "chr9:1-138394717" "chr10:1-133797422" "chr11:1-135086622" "chr12:1-133275309" "chr13:1-114364328" "chr14:1-107043718" "chr15:1-101991189" "chr16:1-90338345" "chr17:1-83257441" "chr18:1-80373285" "chr19:1-58617616" "chr20:1-64444167" "chr21:1-46709983" "chr22:1-50818468")

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
    --compress_fastq\
    --gpu_runners_per_device 15
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Merge fastq files - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' zcat base-calling/fastq${fastq_v}/pass/fastq_runid*.fastq.gz > $fastq
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
  conda activate ngmlr_env
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Alignment step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' ngmlr -t $CPU -r $ref_fasta -q $fastq -x ont \
    |samtools view -Shu |samtools sort -@ $CPU -o $bam --output-fmt BAM
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Bam indexing - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' samtools index $bam -@ $CPU
  conda deactivate
  kill $measure_pid
fi

if [ -e "${vcf}" ]; then
  echo "${vcf} already exists."
else
  measure_resource_usage "${home}" "${inter_med}" "setup${setup_v}-variant_calling.txt" &
  measure_pid=$!
  echo "PID measure $measure_pid" >> "setup${setup_v}-variant_calling.txt" # Escribir el PID en el LOGFILE
  sleep "$((inter_med + 2))"

  echo "----------------"
  echo "Variant calling step"
  conda activate nanopolish
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Fastq indexing step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nanopolish index -d ./raw_reads $fastq
  for region in "${contigs[@]}"; do
    contig=$(echo "$region" | cut -d':' -f1)
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Variant calling step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nanopolish variants \
      --outfile="${chr}/${contig}.vcf" --threads=$CPU\
      --ploidy=2\
      --reads=$fastq\
      --bam=$bam\
      --genome=$ref_fasta\
      --window="${region}"
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compressing vcf files step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bgzip -k -f "${chr}/${contig}.vcf"
    /shared/home/sorozcoarias/anaconda3/bin/time -f 'Indexing vcf files step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' tabix -p vcf "${chr}/${contig}.vcf.gz"
  done
  cd $chr
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Merging vcf files step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bcftools merge chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz --force-samples -o $vcf
  cd $home
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Compressing merged vcf step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' bgzip -k -f "${vcf}"
  /shared/home/sorozcoarias/anaconda3/bin/time -f 'Indexing merged vcf step - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' tabix -p vcf "${vcf}.gz"
  conda deactivate
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

