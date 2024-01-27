patterrnlist=("-1" "0" "1" "2" "3" "4" "5")
issversion=("published2016" "PhyRep2021" "pass7.8")

mkdir -p HighTemplateFit_CCSys
cd HighTemplateFit_CCSys

for issused in $issversion
do 
    for pattern in $patterrnlist
    do
        mkdir -p $issused/p$pattern
        cd $issused/p$pattern
        run_parallel_generic_job --partition c18m --arguments "--binningversion 525version --Rigidityofdataset data_negative --FitMethod ccfixed --issversion $issused --pattern $pattern --ifVGGNN No" --timelimit 5:00 Antiproton_TF2D 0 0
        cd ../..
    done
    mkdir -p $issused/p0VGG
    cd $issused/p0VGG
    run_parallel_generic_job --partition c18m --arguments "--binningversion 525version --Rigidityofdataset data_negative --FitMethod ccfixed --issversion $issused --pattern 0 --ifVGGNN Yes" --timelimit 5:00 Antiproton_TF2D 0 0
    cd ../..
done
cd ..

