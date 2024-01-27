#!/bin/sh

mydir=$1

cp out*.root $mydir/OUTROOT/.
cp out*.log  $mydir/OUTLOG/.
cat finish.dat >> $mydir/finish.dat
cat bad.dat >> $mydir/bad.dat
exit

