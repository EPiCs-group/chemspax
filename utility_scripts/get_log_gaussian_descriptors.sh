#!/bin/bash
# script used to iterate over folders and get some descriptors from gaussian log files
# the first and last SCF energy (used to compare quality of initial geometry)
# the number of optimizations gaussian took
# the total CPU time that the job took in seconds
# HOMO-LUMO gap in hartree (1 hartree = 27.2113961 ev)

FILENAME="bp86.log" # filename of gaussian log file to get values from
N=`wc -l dir| cut -d' ' -f 1`
echo file_name,energy_first,energy_last,n_optimizations,cpu_time,homo_lumo_gap_alpha,homo_lumo_gap_beta
for i in `seq 1 ${N}`; do d=`sed "${i}q;d" dir`
    cd $d
    N_OPTIMIZATIONS=`grep -c "SCF Done" $FILENAME`
    FIRST_TRAJECTORY_ENERGY=`grep "SCF Done" $FILENAME | head -1 | awk '{print $5}'`
    LAST_TRAJECTORY_ENERGY=`grep "SCF Done" $FILENAME | tail -1 | awk '{print $5}'`
    DAYS_CPU_TIME=`grep "Job cpu time" $FILENAME | head -1 | awk '{print $4}'`
    HOURS_CPU_TIME=`grep "Job cpu time" $FILENAME | head -1 | awk '{print $6}'`
    MINUTES_CPU_TIME=`grep "Job cpu time" $FILENAME | head -1 | awk '{print $8}'`
    SECONDS_CPU_TIME=`grep "Job cpu time" $FILENAME | head -1 | awk '{print $10}'`
    TOTAL_CPU_TIME=`convert_gaussian_job_cpu_time.py $DAYS_CPU_TIME $HOURS_CPU_TIME $MINUTES_CPU_TIME $SECONDS_CPU_TIME`
    HOMO_LUMO=`get_homo_lumo_gaussian_log.sh $FILENAME`
    echo $d,$FIRST_TRAJECTORY_ENERGY,$LAST_TRAJECTORY_ENERGY,$N_OPTIMIZATIONS,$TOTAL_CPU_TIME,$HOMO_LUMO
    cd ..
done
