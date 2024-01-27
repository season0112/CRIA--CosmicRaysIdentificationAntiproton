#!/bin/sh

# requires existence of anadir
# inside anadir there are starter.sh
# which does things and its input files
# and ender.sh which handles the outputs
# and requires in input the directory where to put them

# temporary LSF script name 
TMPSCRIPT="job$$_$RANDOM"
SOURCEFILE=/storage/gpfs_ams/ams/users/fnozzoli/amsvar.sh

# is the script existing?
while [ -e "$TMPSCRIPT.slrum" ];
do
echo "Temporary script $TMPSCRIPT.slrum already exists, loop!"
TMPSCRIPT="job$$_$RANDOM"
done    

mv anadir anadir_$TMPSCRIPT

# create the tmp script
cat > $TMPSCRIPT.slrum << EOF
#!/bin/sh
# Automatic script file generated by $0 
#

source $SOURCEFILE
startdir=$PWD
mydir=fnozzoli_$$_$RANDOM
cd /data/
mkdir \$mydir
cd \$mydir
# the next 3 seem to be comments but they are not
#srun $TMPSCRIPT.slrum
#srun -o $PWD/log/$TMPSCRIPT.out
#srun -e $PWD/log/$TMPSCRIPT.err
cp \$startdir/anadir_$TMPSCRIPT/* .
rm -r \$startdir/anadir_$TMPSCRIPT
RUN=\`squeue | grep "RUNNING" | awk '{ print $1 }'\`; echo "the number of runs are \$RUN"
/usr/bin/time -f "%E time-past,\t%S system,\t%U user,\t%P efficiency" ./starter.sh
RUN=\`squeue | grep "run;" | awk '{ print $1 }'\`; echo "the number of runs are \$RUN"
./ender.sh \$startdir
cd ..
rm -Rf \$mydir
exit
EOF

# permissions
chmod 777 $TMPSCRIPT.slrum

# execution of the tmp script 
JOBNAME=`bsub < $TMPSCRIPT.slrum`

# printout
echo "LSF script $TMPSCRIPT.slrum as $JOBNAME"

mv $TMPSCRIPT.slrum $PWD/log/


sleep 3s
