#!/usr/bin/env bash
# redirect stdout/stderr to a file
#exec &> logfile.txt

echo ""
echo "---------------------------------------------------------------------------------------------"
echo "Looks like you want to add some substituents to a metal-ligand complex & optimize 'em."
echo "---------------------------------------------------------------------------------------------"

#rm ../substituents_xyz/automatically_generated/*.xyz
#rm ../substituents_xyz/visualizations/*.png
#
CURRENT_DIR=$(pwd)
SKELETON_LIST=$(cd skeletons && ls | cut -d '.' -f 1)
# select 1 random substituent with C as central atom as starting point
STARTING_C_SUBSTITUENT=$(cd substituents_xyz/manually_generated/ && ls -d C* | xargs shuf -n1 -e | cut -d '.' -f 1)
# select 5 random substituents with C as central atom
RANDOM_C_SUBSTITUENTS=$(cd substituents_xyz/manually_generated/ && ls -d C* | xargs shuf -n5 -e | cut -d '.' -f 1)

#SOURCE_FILE="/skeletons/RuPNP_iPr_skl"
#SOURCE_NAME=$(echo ${SOURCE_FILE} | cut -d '/' -f 3)
#TARGET_NAME=${SOURCE_NAME}_func

for skeleton in ${SKELETON_LIST}; do
# loop over skeleton list and set index back to 1
    SOURCE_FILE=skeletons/${skeleton}
    TARGET_NAME=${skeleton}_func
    i=1
    echo "creating initial file"
    # functionalize and optimize initial functionalized version of skeleton
    python3 main_attach_substituent.py ${SOURCE_FILE}.xyz ${TARGET_NAME}_${i} ${STARTING_C_SUBSTITUENT} substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.54
    # optimization
    cd substituents_xyz/automatically_generated/
    xtb ${TARGET_NAME}_${i}.xyz --opt --chrg 0 --uhf 0 --gbsa acetonitrile > xtb.out
    # clean up mess and move relevant file to correct folder
    mv xtbopt.xyz optimized_structures/${TARGET_NAME}_${i}_opt.xyz
    rm -f xtbrestart
    cd -
        for sub in ${RANDOM_C_SUBSTITUENTS}; do
        echo "Running recursive loop, run:" ${i}
        python3 main_attach_substituent.py substituents_xyz/automatically_generated/${TARGET_NAME}_${i}.xyz ${TARGET_NAME}_$((i+1)) ${sub} substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.54
        # optimization
	cd substituents_xyz/automatically_generated
        xtb ${TARGET_NAME}_$((i+1)).xyz --opt --chrg 0 --uhf 0 --gbsa acetonitrile > xtb.out
        # clean up mess and move relevant file to correct folder
        mv xtbopt.xyz optimized_structures/${TARGET_NAME}_$((i+1))_opt.xyz
        rm -f xtbrestart
	cd -
        i=$((i+1))
        done
done
