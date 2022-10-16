# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
import glob
import subprocess
""" Run script to convert all .gjf files in current path to .xyz files
"""


def remove_last_line(filename):
    with open(filename) as f:
        lines = f.readlines()
        last = len(lines) - 1
        lines[last] = lines[last].replace('\r', '').replace('\n', '')
    with open(filename, 'w') as wr:
        wr.writelines(lines)


def count_atoms_and_write(filename):
    lines = open(filename).readlines()
    count = len(lines)
    count -= 2
    lines[0] = str(count)+'\n'
    with open(filename, 'w') as wr:
        wr.writelines(lines)


def convert_gjf_to_xyz_processing():
    for file in glob.glob("*.xyz"):
        remove_last_line(file)
        count_atoms_and_write(file)


if __name__ == "__main__":
    subprocess.call("gjf_to_xyz.sh", shell=True)
    convert_gjf_to_xyz_processing()
