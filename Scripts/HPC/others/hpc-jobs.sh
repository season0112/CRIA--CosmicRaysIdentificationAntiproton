#!/usr/bin/env bash

##### Old project:jara0052, New project:p0020149
while [ 1 ]
do
    echo Sichen: $(squeue| grep 'bo791269'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'bo791269'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Robin: $(squeue| grep 'rs429310'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'rs429310'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    #echo Fabian: $(squeue| grep 'fm606351'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'fm606351'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Yasaman: $(squeue| grep 'op115134'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'op115134'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Emiliano: $(squeue| grep 'ms255555'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'ms255555'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo ChanHoon: $(squeue| grep 'cc435041'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'cc435041'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Manbing: $(squeue| grep 'lk974211'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'lk974211'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Sophie: $(squeue| grep 'hh288529'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'hh288529'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\) 
    #echo Nicolay: $(squeue| grep 'nn058590'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'nn058590'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Henning: $(squeue| grep 'hg781119'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'hg781119'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Valery: $(squeue| grep 'vz897724'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'vz897724'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    #echo Leila: $(squeue| grep 'la459195'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'la459195'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo Schael: $(squeue| grep 'ss243125'|grep 'RUNNING'|wc|awk '{print $1}')/$(squeue| grep 'ss243125'|grep 'PENDING'|wc|awk '{print $1}') \(RUNNING/PENDING\)
    echo ## echo means skip a line. 
    sleep 3
done


#echo $USER: $(bjobs -r -noheader 2>/dev/null | wc -l)/$(bjobs -noheader 2>/dev/null | wc -l)
#echo AMS: $(bjobs -r -u all -P jara0052 -noheader 2>/dev/null | wc -l)/$(bjobs -u all -P jara0052 -noheader 2>/dev/null | wc -l)
#echo JARA: $(bjobs -r -u all -q jara-clx -noheader 2>/dev/null | wc -l)/$(bjobs -u all -q jara-clx -noheader 2>/dev/null | wc -l)
#echo ALL: $(bjobs -r -u all -noheader 2>/dev/null | wc -l)/$(bjobs -u all -noheader 2>/dev/null | wc -l)
#echo GPU: $(bhosts -l claix-gpu | grep closed | wc -l)/$(lshosts | grep lng | wc -l)


