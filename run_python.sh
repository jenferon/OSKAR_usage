#! /usr/bin/bash
#SBATCH -p shortq,defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=1:00:00
#SBATCH --output=log/%x.%j.o
#SBATCH --error=log/%x.%j.e
#SBATCH --array=0

source ~/.bash_profile

cd /gpfs01/home/ppxjf3/OSKAR/
activate jen

python testing_noise.py
