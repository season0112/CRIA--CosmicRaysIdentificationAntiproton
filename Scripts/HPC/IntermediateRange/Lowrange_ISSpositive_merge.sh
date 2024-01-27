#cd $1

#for i in $(seq 0 $2)
#do
#    var=$(printf "%05d" "$i")
#    echo Now it is submitting: $var;
#    run_parallel_generic_job.py --keep_directories --arguments "--input AntiprotonIntermediateEnergy_Tree_${var}_* --resultprefix reduced_root/AntiprotonIntermediateEnergy_Tree_${var}" ac_merge 0 0;
#done

#run_parallel_generic_job.py --arguments "--input AntiprotonIntermediateEnergy_Tree_%05d_* --resultprefix reduced_root/AntiprotonIntermediateEnergy_Tree_%05d" ac_merge 0 $2
run_parallel_generic_job.py --arguments "--input $1/AntiprotonIntermediateEnergy_Tree_%05d_* --resultdir $1/reduced_root/ --resultprefix AntiprotonIntermediateEnergy_Tree_%05d" ac_merge 0 999

#while [ $(squeue | grep 'bo791269' | grep 'ac_merge'|wc |awk '{print $1}') != 0 ]
#do
#    echo some jobs unfinished, sleeping for 2 minutes...
#    sleep 2m
#done 


#find . -name "AntiprotonIntermediateEnergy_Tree*"| cut -c 37-41|sort> ../../reduced_index.txt
#cd ../../
#for i in $(seq 0 999); do printf "%05d\n" "$i"; done > reduced_all.txt
#diff reduced_all.txt reduced_index.txt| awk '/</' | cut -c 3-7> reduce_missing_index.txt


#for j in $(cat reduce_missing_index.txt)
#do 
#    cd $1
#    run_parallel_generic_job.py --keep_directories --arguments "--input AntiprotonIntermediateEnergy_Tree_${j}_* --resultprefix reduced_root/AntiprotonIntermediateEnergy_Tree_${j}" ac_merge 0 0;
#    cd ..
#done

