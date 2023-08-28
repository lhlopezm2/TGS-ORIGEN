#!/usr/bin/bash
#SBATCH -J nanocall
#SBATCH -D .
#SBATCH -o results/out-nanocall.txt
#SBATCH -e results/err-nanocall.txt
#SBATCH -n 22
#SBATCH -N 1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:3g.20gb:1
#SBATCH --account=coffea_genomes
#SBATCH --time=1-23:00:00
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luish.lopezm@autonoma.edu.co
source ~/.bashrc
CPU=20
conda activate nanocall
/shared/home/sorozcoarias/anaconda3/bin/time -f 'Bam indexing - Elapsed Time: %e s - Memory used: %M kB -CPU used: %P' nanocall ./raw_reads.fofn -t $CPU >output.fa
conda deactivate
