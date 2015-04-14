#!/bin/bash

# Set coordinates type coor (1 for Mercator, 2 for lat-lon)
coor=2
sed -i "s/^coor *= *[0-9]\+/coor=$coor/" control_file.txt

# Run igw_v1 for all following vals of nph
for nph in 3600 4320 #5040
do
sed -i "s/^nph *= *[0-9]\+/nph=$nph/" control_file.txt
unbuffer ./igw_v1 2>&1 | tee -a logs/itd_runs/log_${coor}_${nph}.txt
done

