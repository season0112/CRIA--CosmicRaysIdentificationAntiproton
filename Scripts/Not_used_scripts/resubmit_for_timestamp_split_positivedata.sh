cd $1/output
grep -nLr 'exited with return code 0'|sort |cut -c 8-12>../jobs/missing_index.txt
cd ../jobs/

for i in $(cat missing_index.txt)
do
    var=$(printf "%05d" "$i")
    sbatch slurmjob_${var}.sh
done    

