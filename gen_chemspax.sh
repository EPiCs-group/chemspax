#!/bin/bash
echo "correct usage is sh gen_chemspax.sh skeleton_name number_of_substituent_sites"

mkdir "slurmfiles_$1"

N=`wc -l substitutes_list/substitutes_combined_$2.txt| cut -d' ' -f 1`

for j in `seq 1 $N`;do
cat > slurmfiles_$1/chemspax_${1}_$j.slurm <<EOF
#!/bin/bash
#SBATCH --job-name=${1}${j}chemspax
#SBATCH --mem-per-cpu=500MB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --time=04:59:59
EOF
d=`sed "${j}q;d" substitutes_list/substitutes_combined_$2.txt`
filename=$1_$d
filename="${filename// /_}"
filename="${filename//\"/}"
echo "sh functionalize_and_optimize_obabel.sh C $d" >> slurmfiles_$1/chemspax_${1}_$j.slurm
sbatch slurmfiles_$1/chemspax_${1}_${j}.slurm
done