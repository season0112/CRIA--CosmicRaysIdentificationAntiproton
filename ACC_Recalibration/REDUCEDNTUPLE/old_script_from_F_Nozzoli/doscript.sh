#!/bin/sh
rm -r log/*
rm -r OUTLOG/*
rm -r anadi*
rm -r OUTROOT/*
rm bad.dat
touch bad.dat
rm finish.dat
touch finish.dat

totrun=$(cat listarun.dat | wc -l)
fatrun=0
dpat=$(cat path.pat)
dout="out"
while read line
do
djob=$(bjobs | grep PEND | wc -l)
while [ $djob -gt 200 ]
do
sleep 60s
djob=$(bjobs | grep PEND | wc -l)
done

rm -r anadir
mkdir anadir
cp *.C anadir/.
cp *.h anadir/.
cp Makefile anadir/.
touch anadir/finish.dat
touch anadir/bad.dat
cp doscript_ender.sh anadir/ender.sh
chmod 777 anadir/ender.sh
cat > anadir/starter.sh << EOF
#!/bin/sh
# Automatic script file generated by $0
#
source /storage/gpfs_ams/ams/users/fnozzoli/amsvar.sh
make
./analyze $dpat$line $dout$line
exit
EOF

chmod 777 anadir/starter.sh
./dojob_new.sh
# anadir rinominata e poi cancellata in dojob_new
fatrun=$(($fatrun+1))
echo "submitted " $fatrun " of " $totrun 
done < listarun.dat
exit

