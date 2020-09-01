#!/usr/bin/env bash
# redirect stdout/stderr to a file
exec &> logfile.txt

echo ""
echo "-------------------------------------------------------"
echo "Looks like you want to add some substituents to a metal-ligand complex."
echo "-------------------------------------------------------"


SOURCE_FILE="skeletons/RuPNP_iPr_skl"
TARGET_FILE="substituents_xyz/automatically_generated/RuPNP_CH4"

rm substituents_xyz/automatically_generated/*.xyz
rm visualizations/*.png

echo "creating initial file"
python main.py ${SOURCE_FILE}.xyz RuPNP_CH4_1 None recursive H H H False False

for i in `seq 1 3`
do
  echo "Running recursive loop, #run:" ${i}
  python main.py ${TARGET_FILE}_${i}.xyz RuPNP_CH4_$((i+1)) None recursive H H H False False
done