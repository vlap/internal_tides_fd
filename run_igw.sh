#!/bin/bash

USER="amtvl"
MODEL="internal_fd"

if [ "$HOSTNAME" == "theon" ]; then
   DST="/scratch/1/users/"$USER"/Tides/data/LAG/"$MODEL"/out/"
   DRBOX="/scratch/1/users/"$USER"/Dropbox/work/runs/"$MODEL"/"
elif [ "$HOSTNAME" == "xigor" ]; then
   DST="/scratch/a/users/"$USER"/Tides/data/LAG/"$MODEL"/out/"    
else
   echo "Modify the script. Only works of theon and xigor"
   exit 1
fi

 # assume that there is no significant delay in generation of the time-stamp
   NOW=$(date +"%Y_%m_%d__%H_%M")

nice -n 19 unbuffer ./igw_v1 2>&1 | tee -a $DST/log$NOW.txt # $DRBOX/$NOW/log.txt
mv $DST/log$NOW.txt $DST/$NOW/

# mkdir $DST/$NOW
