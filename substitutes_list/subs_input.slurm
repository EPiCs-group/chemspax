#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=180GB
#SBATCH --partition=
#SBATCH --time 23:59:59
#SBATCH -J substituents

python substitutes_combined.py
