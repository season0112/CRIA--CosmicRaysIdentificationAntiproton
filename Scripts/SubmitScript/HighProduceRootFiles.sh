patterrnlist=("-1" "0" "1" "2" "3" "4" "5")

#issversion=("published2016" "PhyRep2021" "pass7.8")
issversion=("published2016" "PhyRep2021")

mkdir -p HighProduceRootFilesSubmit
cd HighProduceRootFilesSubmit

for issused in $issversion
do
    for pattern in $patterrnlist
    do
        mkdir -p $issused/p$pattern
        cd $issused/p$pattern
        run_parallel_generic_job --partition c18m --arguments "--issversion $issused --binningversion 525version --pattern $pattern --ifVGGNN No" --timelimit 5:00 $MY_ANALYSIS/Scripts/HPC/HighRange/antip_p_Number_TH2D.py 0 0
        cd ../..
    done
    mkdir -p $issused/p0VGG
    cd $issused/p0VGG
    run_parallel_generic_job --partition c18m --arguments "--issversion $issused --binningversion 525version --pattern 0 --ifVGGNN Yes" --timelimit 5:00 $MY_ANALYSIS/Scripts/HPC/HighRange/antip_p_Number_TH2D.py 0 0
    cd ../..
done
cd ..




