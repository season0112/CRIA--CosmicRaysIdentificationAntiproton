#!/bin/bash

cd $JUAMSHIGHENERGYDATADIR/$1/$2/results
mkdir -p lists
find $PWD -name "*_Tree_*" > lists.list
mv lists.list lists
cd lists
split -d -a 5 -l $3 lists.list x
mkdir -p filelists
mv x* filelists

