patterrnlist=("0" "1" "2" "4")
issversion=("published2016" "PhyRep2021" "pass7.8")
RigiditySign=("negative" "positive")

mkdir -p HighTemplateFit_forAllPatternsSubmit
cd HighTemplateFit_forAllPatternsSubmit

for sign in $RigiditySign
do 
    mkdir -p $sign
    cd $sign
    for issused in $issversion
    do
        mkdir -p $issused
        cd $issused
        for pattern in $patterrnlist
        do
            mkdir -p p$pattern
            cd p$pattern
            run_parallel_generic_job --partition c18m --arguments "--binningversion 525version --Rigidityofdataset data_$sign --FitMethod ccfree --issversion $issused --pattern $pattern --ifVGGNN No" --timelimit 5:00 Antiproton_TF2D 0 0
            cd ..
        done
        mkdir -p p0VGG
        cd p0VGG
        run_parallel_generic_job --partition c18m --arguments "--binningversion 525version --Rigidityofdataset data_$sign --FitMethod ccfree --issversion $issused --pattern 0 --ifVGGNN Yes" --timelimit 5:00 Antiproton_TF2D 0 0
        cd .. 
        cd ..
    done
    cd ..
done

cd ..


