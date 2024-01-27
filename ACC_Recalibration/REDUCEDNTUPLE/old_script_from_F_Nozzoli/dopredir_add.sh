#!/bin/sh


echo "Run in area update for path data analysis:"
cat path.pat
echo " " 

dpat=$(cat path.pat)
rm listappo.dat
rm listarun.dat
rm listabadrun.dat
ls $dpat | grep root >> listappo.dat 

while read line
do
tst=$(cat $AMSDataDir/BadRuns/TRD_U*:* | grep ${line:0:10} | wc -l)
tst2=$(cat $AMSDataDir/BadRuns/TRD_B*:* | grep ${line:0:10} | wc -l)
if [ $tst -eq 0 ]; then
if [ $tst2 -eq 0 ]; then
echo $line >> listarun.dat
fi
fi
if [ ! $tst -eq 0 ]; then
echo $line >> listabadrun.dat
fi
if [ ! $tst2 -eq 0 ]; then
echo $line >> listabadrun.dat
fi
done < listappo.dat

touch finish.dat
touch bad.dat

exit

