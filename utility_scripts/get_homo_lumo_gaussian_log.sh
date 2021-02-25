#!/bin/bash
# script to get the homo and lumo from gaussian log files
# made by Vivek Sinha
H=`awk '/Alpha  occ. eigenvalues/ {print $NF}' $1 | tail -1`
L=`awk '/Alpha virt. eigenvalues/ {print $5}' $1 | head -1`
if grep -q "Beta  occ. eigenvalues" $1
        then
        Hb=`awk '/Beta  occ. eigenvalues/ {print $NF}' $1 | tail -1`
        Lb=`awk '/Beta virt. eigenvalues/ {print $5}' $1 | head -1`
        echo "Alpha:" $H,$L, "Beta:" $Hb,$Lb
else
        echo $H,$L
fi
