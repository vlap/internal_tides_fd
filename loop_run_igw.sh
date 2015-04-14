#!/bin/bash

# Set coordinates type coor (1 for Mercator, 2 for lat-lon)
coor=2
sed -i "s/^coor *= *[0-9]\+/coor=$coor/" control_file.txt


# Run igw_v1 for all following vals of nph
for nph in 3600 4500
do
sed -i "s/^nph *= *[0-9]\+/nph=$nph/" control_file.txt
./run_igw.sh # | tee -a logs/itd_runs/log_${coor}_${nph}.txt
done

#n_modes=3
#sed -i "s/^n_modes *= *[0-9]\+/n_modes=$n_modes/" control_file.txt

#N_hmin = 50
#sed -i "s/^N_hmin *= *[0-9]\+/N_hmin=$N_hmin/" control_file.txt

#for smooth_type in 1 3
#do
#sed -i "s/^smooth_type *= *[0-9]\+/smooth_type=$smooth_type/" control_file.txt
#unbuffer ./igw_v1 2>&1 | tee -a logs/log_${coor}_${nph}_${smooth_type}.txt
#done

#nph=3600
#for itm_qmax in 1d-3 1d-2 1d-1 1d+0
#do
#sed -i "s/^itm_qmax *= *[0-9]\+d.[0-9]\+/itm_qmax=$itm_qmax/" control_file.txt
#unbuffer ./igw_v1 2>&1 | tee -a logs/log_${coor}_${nph}_${smooth_type}.txt
#done
