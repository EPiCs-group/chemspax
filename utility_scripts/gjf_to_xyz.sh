#!/usr/bin/env bash
# https://askubuntu.com/questions/1184243/how-to-convert-gjf-format-to-xyz-format
# run from gjf_to_xyz.py to use this script on all .gjf files in current path
exec &> logfile.txt # redirect stdout/stderr to a file

for file_name in *.gjf; do
tail -1 ${file_name}  > ${file_name%.*}.xyz
echo"" >> ${file_name%.*}.xyz
grep '[0-9]\.[0-9][0-9]' ${file_name} >> ${file_name%.*}.xyz
done