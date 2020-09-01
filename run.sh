#!/usr/bin/env bash
# redirect stdout/stderr to a file
exec &> logfile.txt

echo ""
echo "-------------------------------------------------------"
echo "Looks like you want to add some substituents to a metal-ligand complex."
echo "-------------------------------------------------------"


SOURCE_FILE="skeletons/RuPNP_iPr.xyz"
TARGET_FILE="substituents_xyz/automatically_generated/RuPNP_CH4.xyz"
i=1

if [ -f ${TARGET_FILE} ]
then
    echo ${TARGET_FILE} "exists already. Removing it first."
    rm -rf ${TARGET_FILE}
    cp skeletons/RuPNP_iPr.xyz ${TARGET_FILE}
else
    cp skeletons/RuPNP_iPr.xyz ${TARGET_FILE}
fi

echo "creating initial file"
python main.py ${SOURCE_FILE} RuPNP_CH4 [[19,4],[23,4],[24,3],[20,3]] recursive H H C False False

until [ $i -gt 3 ]
do
  echo "Running recursive loop, #run:" ${i}
  python main.py ${TARGET_FILE} RuPNP_CH4 None recursive H H C False False
  ((i=i+1))
done


#else
#    echo "Running loop"
#    python main.py substituents_xyz/manually_generated/CH3.xyz CH4 [[1,0],[4,0]] initial H H C False False
#    python main.py substituents_xyz/automatically_generated/CH4.xyz CH4 None recursive H H H True False
#fi