#!/bin/bash

#while read line
#do
#  echo $line | xargs echo -n
#  echo $line
#  echo $line |xargs echo -n
#  cd OUTROOT_nominal
#  tem=$($line |xargs echo -n)
#  echo $tem
#  cp $tem ../gap_rootfile
#  cd ..

#done < gap.txt


for i in {1..5978}
do 
  tem=$(sed -n $((i)),$((i))p gap2.dat)
  cd OUTROOT_nominal
  cp $tem ../gap_rootfile2
  cd ..
done



