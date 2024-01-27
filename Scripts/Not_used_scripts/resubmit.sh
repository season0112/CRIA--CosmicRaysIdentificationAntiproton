#!/usr/bin/env bash

#### submit unstarted jobs:
echo submit unstarted jobs now.
cd $1/$2/$3/output/
outputnumber=$(ls|wc|awk '{print $1}');
cd ../jobs/
jobsnumber=$(ls|wc|awk '{print $1}');
if [ $aa != $bb ]; then
    cd $1/$2/$3/output
    ls | cut -c8-12 > ../outtest.txt
    cd ../jobs
    rm core.*
    ls | cut -c10-14 > ../jobtest.txt
    cd ..
    diff outtest.txt jobtest.txt| awk '/>/' | cut -c3-7>indexdiff.txt
    for number in $(cat indexdiff.txt);
      do
        cd jobs
        sbatch slurmjob_${number}.sh
        cd ..
    done
fi
echo all jobs started.

echo checking if all jobs finised...
cd $1/$2/$3/output/
aaa=$(grep -nr 'exited with return code'|wc | awk '{print $1}') #all complete jobs.
bbb=$(grep -nr 'CANCELLED AT'|wc | awk '{print $1}') ## canceled jobs.
cd ../jobs
ccc=$(ls|wc|awk '{print $1}')
while [ $(($aaa+$bbb)) != $ccc ] 
do 
    echo some unfinished jobs, wait for 2 mins...
    sleep 2m
done
echo all jobs finised.

#### submit crashed jobs:
echo checking if there are some crashed jobs...
cd $1/$2/$3/output
echo $(grep -nrL 'exited with return code 0'|wc)> ../temporary.txt
while [ $(awk '{print $1}' ../temporary.txt) != 0 ] 
do
    echo there are some crashed jobs. resubmit now.
    cd $1
    cd $2/$3/output/
    print $PWD
    grep -nrL 'exited with return code 0' > $2$3_tem.txt
    cut -c8-12 ./$2$3_tem.txt > $2$3_tem2.txt
    #grep -nr 'exited with return code 134' > $2$3_tem.txt
    #awk '{print $2}' ./$2$3_tem.txt > $2$3_tem2.txt
    rm $2$3_tem.txt
    mv $2$3_tem2.txt ../../../
    cd ../../..
    #cd ..
    print $PWD
    for number in $(cat $2$3_tem2.txt);
      do 
        cd $2/$3/jobs/;
        print slurmjob_${number}.sh;
        sbatch slurmjob_${number}.sh;
        cd ../../..;
    done   
    rm $2$3_tem2.txt

    cd $1/$2/$3/output/
    while [ $(ls|wc|awk '{print $1}') != $(echo $(grep -nr 'exited with return'|wc) | awk '{print $1}') ]
    do
        echo waiting for resubmitted jobs...
        sleep 2m
    done
done
echo all resubmitted jobs finished. 

cd $1/$2/$3/
rm temporary.txt
cd $MY_ANALYSIS/Scripts/HPC/others/

