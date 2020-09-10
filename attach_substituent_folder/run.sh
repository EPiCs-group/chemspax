#!/usr/bin/env bash
# redirect stdout/stderr to a file
exec &> logfile.txt

echo ""
echo "-------------------------------------------------------"
echo "Looks like you want to add some substituents to a metal-ligand complex."
echo "-------------------------------------------------------"


SOURCE_FILE="../skeletons/RuPNP_iPr_skl"
TARGET_NAME="RuPNP_iPr_sub"
SUBSTITUENT_LIST="CCl3 CH2F CF3 CH2Cl CH2OCH3"
rm ../substituents_xyz/automatically_generated/*.xyz
rm ../substituents_xyz/visualizations/*.png

echo "creating initial file"
python ../main_attach_substituent.py ${SOURCE_FILE}.xyz ${TARGET_NAME}_1 CH3 ../substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.54

i=1
for sub in ${SUBSTITUENT_LIST}; do
    echo "Running recursive loop, run:" ${i}
    python ../main_attach_substituent.py ../substituents_xyz/automatically_generated/${TARGET_NAME}_${i}.xyz ${TARGET_NAME}_$((i+1)) ${sub} ../substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.54
    i=$((i+1))
done
