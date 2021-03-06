#!/bin/bash

USER="amtvl"
MODEL="internal_fd"

if [ "$HOSTNAME" == "theon" ]; then
   HOST="xigor"
   DATADIRSRC="/scratch/a/users/"$USER"/Tides/data/LAG/"$MODEL"/out/"
   DST="/scratch/1/users/"$USER"/Tides/data/LAG/"$MODEL"/out/m2_only/"
elif [ "$HOSTNAME" == "xigor" ]; then
   HOST="theon"
   DATADIRSRC="/scratch/1/users/"$USER"/Tides/data/LAG/"$MODEL"/out/"
   DST="/scratch/a/users/"$USER"/Tides/data/LAG/"$MODEL"/out/"   
else
   echo "Modify the script. Only works of theon and xigor"
   exit 1
fi

EXCLUDE1="temp"
#EXCLUDE2="mats"

SSH_COMMAND="ssh "
EXTRA="$2" #"--dry-run" "-L" # copy files that links refer to

#Run directory name
DIR="$1"
if [ -n "$DIR" ]; then

    if [ "$DIR" == "last" ]; then 
        SRC=$USER@$HOST:$DATADIRSRC"$(ssh $HOST 'ls -tr '$DATADIRSRC' | tail -1')" # copy the last run
    else
        SRC=$USER@$HOST:$DATADIRSRC$DIR  # Name of the run to be copied
    fi
else 
    exit 1
fi

# Remote backup. First establish ssh connection to amsta. Set permission rwX to all
rsync -avh $EXTRA --progress --exclude "$EXCLUDE1" --exclude "$EXCLUDE2" -e "$SSH_COMMAND" $SRC $DST # --stats --human-readable

# Some of the symbolic links will be broken now (they have absolute paths). Adjust the paths
find $DST -lname $DATADIRSRC'*' -print | while read file
do
    TARGET=$(readlink $file)
    NEW=$(echo $TARGET | sed 's:^'$DATADIRSRC':'$DST':')
    rm $file
    ln -s $NEW $file
done
