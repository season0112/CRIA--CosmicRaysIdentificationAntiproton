cd $HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/

mkdir -p binmerge1;
mkdir -p binmerge2;

mv *2.97_3.29* binmerge1;
mv *3.29_3.64* binmerge1;
mv *3.64_4.02* binmerge1;
mv *4.02_4.43* binmerge1;
mv *4.43_4.88* binmerge1;
mv *4.88_5.37* binmerge1;
mv *5.37_5.9* binmerge1;
mv *5.9_6.47* binmerge1;
mv *6.47_7.09* binmerge1;
mv *7.09_7.76* binmerge1;
mv *7.76_8.48* binmerge1;
mv *8.48_9.26* binmerge1;
mv *9.26_10.1* binmerge1;
mv *10.1_11* binmerge1;
mv *11_12* binmerge1;
mv *12_13* binmerge1;
mv *13_14.1* binmerge1;
mv *14.1_15.3* binmerge1;
mv *15.3_16.6* binmerge1;
mv *16.6_18* binmerge1;


mv *2.97_3.64* binmerge2;
mv *3.64_4.43* binmerge2;
mv *4.43_5.37* binmerge2;
mv *5.37_6.47* binmerge2;
mv *6.47_7.76* binmerge2;
mv *7.76_9.26* binmerge2;
mv *9.26_11* binmerge2;
mv *11_13* binmerge2;
mv *13_15.3* binmerge2;
mv *15.3_18* binmerge2;

mv *binmerge1.* binmerge1
mv *binmerge2.* binmerge2


