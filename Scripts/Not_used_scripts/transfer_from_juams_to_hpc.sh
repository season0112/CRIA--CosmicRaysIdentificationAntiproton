#!/bin/bash
#ssh copy
#cd /hpcwork/jara0052/sichen/analysis_7.0/
#echo $0

#:echo "ls -l; echo 'Hello World'" | ssh copy /bin/zsh

#mkdir myfolder  | ssh copy /bin/zsh

ssh copy "mkdir -p $HPCHIGHENERGYDATADIR/$1/$2"
scp -r $JUAMSHIGHENERGYDATADIR/$1/$2/results/rawdata/transferdata copy:$HPCHIGHENERGYDATADIR/$1/$2/

# scp -r $1/$2/results/rawdata/transferdata copy:/hpcwork/jara0052/sichen/analysis_7.0/$1/$2/





