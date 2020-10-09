#!/usr/bin/env bash
# redirect stdout/stderr to a file
# exec &> logfile.txt

echo ""
echo "---------------------------------------------------------------------------------------------"
echo "Looks like you want to add some substituents to a metal-ligand complex & optimize 'em."
echo "Created by: Adarsh Kalikadien & Vivek Sinha"
echo "---------------------------------------------------------------------------------------------"

functionalize_skeletons_C_substituents (){
SKELETON_LIST=$(cd skeletons && ls | cut -d '.' -f 1)
# select 1 random substituent with C as central atom as starting point
STARTING_C_SUBSTITUENT=$(cd substituents_xyz/manually_generated/ && ls -d C* | xargs shuf -n1 -e | cut -d '.' -f 1)
#STARTING_C_SUBSTITUENT="F"
# select 5 random substituents with C as central atom
RANDOM_C_SUBSTITUENTS=$(cd substituents_xyz/manually_generated/ && ls -d C* | xargs shuf -n5 -e | cut -d '.' -f 1)
#RANDOM_C_SUBSTITUENTS="F"
echo ${STARTING_C_SUBSTITUENT} > substituents_xyz/automatically_generated/substituents.txt
echo ${RANDOM_C_SUBSTITUENTS} >> substituents_xyz/automatically_generated/substituents.txt
# C-C bond length = 1.54 A
# https://phys.org/news/2018-03-carbon-carbon-bond-length.html

# activate option to use backslash options in echo function
shopt -s xpg_echo

for skeleton in ${SKELETON_LIST}; do
# loop over skeleton list and set index back to 1
    SOURCE_FILE=skeletons/${skeleton}
    TARGET_NAME=${skeleton}_func
    i=1
    echo "creating initial file of" ${skeleton} ${STARTING_C_SUBSTITUENT}
    # functionalize and optimize initial functionalized version of skeleton
    python3 main_attach_substituent.py ${SOURCE_FILE}.xyz ${TARGET_NAME}_${i} ${STARTING_C_SUBSTITUENT} substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.54
    # optimization
    cd substituents_xyz/automatically_generated/
    # write functionalization list to optimized file to be able to use that file as source for new functionalizations
    FUNCTIONALIZATION_LIST=$(sed '2q;d' ${TARGET_NAME}_${i}.xyz)
    sed -i '1s/.*/'"${FUNCTIONALIZATION_LIST}"'/' ${TARGET_NAME}_${i}.mol
    cd -
        for sub in ${RANDOM_C_SUBSTITUENTS}; do
        echo "Running recursive loop, run:" ${i} ${sub}
        python3 main_attach_substituent.py substituents_xyz/automatically_generated/${TARGET_NAME}_${i}.xyz ${TARGET_NAME}_$((i+1)) ${sub} substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.54
        # optimization
	    cd substituents_xyz/automatically_generated
        # write functionalization list to optimized file to be able to use that file as source for new functionalizations
        FUNCTIONALIZATION_LIST=$(sed '2q;d' ${TARGET_NAME}_$((i+1)).xyz)
        sed -i '1s/.*/'"${FUNCTIONALIZATION_LIST}"'/' ${TARGET_NAME}_$((i+1)).mol
	    cd -
        i=$((i+1))
        done
done
}

functionalize_skeletons_P_substituents (){
SKELETON_LIST=$(cd skeletons && ls | cut -d '.' -f 1)
# select 1 random substituent with P as central atom as starting point
STARTING_P_SUBSTITUENT=$(cd substituents_xyz/manually_generated/ && ls -d P* | xargs shuf -n1 -e | cut -d '.' -f 1)
# select 5 random substituents with P as central atom
RANDOM_P_SUBSTITUENTS=$(cd substituents_xyz/manually_generated/ && ls -d P* | xargs shuf -n5 -e | cut -d '.' -f 1)
echo ${STARTING_P_SUBSTITUENT} > substituents_xyz/automatically_generated/substituents.txt
echo ${RANDOM_P_SUBSTITUENTS} >> substituents_xyz/automatically_generated/substituents.txt
# C-P bond length = 1.87 A
# https://sites.google.com/site/chempendix/bond-lengths

# activate option to use backslash options in echo function
shopt -s xpg_echo

for skeleton in ${SKELETON_LIST}; do
# loop over skeleton list and set index back to 1
    SOURCE_FILE=skeletons/${skeleton}
    TARGET_NAME=${skeleton}_func
    i=1
    echo "creating initial file of" ${skeleton} ${STARTING_P_SUBSTITUENT}
    # functionalize and optimize initial functionalized version of skeleton
    python3 main_attach_substituent.py ${SOURCE_FILE}.xyz ${TARGET_NAME}_${i} ${STARTING_P_SUBSTITUENT} substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.87
    # optimization
    cd substituents_xyz/automatically_generated/
    # write functionalization list to optimized file to be able to use that file as source for new functionalizations
    FUNCTIONALIZATION_LIST=$(sed '2q;d' ${TARGET_NAME}_${i}.xyz)
    sed -i '1s/.*/'"${FUNCTIONALIZATION_LIST}"'/' ${TARGET_NAME}_${i}.mol
    cd -
        for sub in ${RANDOM_P_SUBSTITUENTS}; do
        echo "Running recursive loop, run:" ${i} ${sub}
        python3 main_attach_substituent.py substituents_xyz/automatically_generated/${TARGET_NAME}_${i}.xyz ${TARGET_NAME}_$((i+1)) ${sub} substituents_xyz/manually_generated/central_atom_centroid_database.csv 1.87
        # optimization
	    cd substituents_xyz/automatically_generated
        # write functionalization list to optimized file to be able to use that file as source for new functionalizations
        FUNCTIONALIZATION_LIST=$(sed '2q;d' ${TARGET_NAME}_$((i+1)).xyz)
        sed -i '1s/.*/'"${FUNCTIONALIZATION_LIST}"'/' ${TARGET_NAME}_$((i+1)).mol
	    cd -
        i=$((i+1))
        done
done
}

# make script runnable with system arguments
case $1 in
    "C")
    functionalize_skeletons_C_substituents
    echo ""
    echo ""
    echo "Successfully functionalized all skeletons with C substituents"

    ;;
    "P")
    functionalize_skeletons_P_substituents
    echo ""
    echo ""
    echo "Successfully functionalized all skeletons with P substituents"

    ;;
    "CP")
    functionalize_skeletons_C_substituents
    functionalize_skeletons_P_substituents
    echo ""
    echo ""
    echo "Sucessfully functionalized all skeletons with P and C substituents"

    ;;
    *)
    echo 'available options: C, P, CP'
    echo 'example usage: bash functionalize_and_optimize_obabel.sh CP'
    exit 1
    ;;
esac
exit 0;
#bash functionalize_skeletons_C_substituents C