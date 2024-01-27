mkdir -p Submit_Antiproton_Low_MakeTimeDependentTemplate
cd Submit_Antiproton_Low_MakeTimeDependentTemplate

for (( i=2; i<18; i=i+2))
do
    echo Now is: rigidity_start = $i ,  rigidity_end = $[i+1]

    left=$i
    right=$[i+1]

    mkdir -p $left"_"$right
    cd $left"_"$right
    
    mkdir -p 6months
    cd 6months
    for (( p=0; p<21; p=p+1))
    do
        LeftTimeIndex=$p
        RightTimeIndex=$[p+1] 
        mkdir -p $LeftTimeIndex"_"$RightTimeIndex
        cd $LeftTimeIndex"_"$RightTimeIndex
        run_parallel_generic_job --partition c18m --timelimit 24:00 --memory 180 --arguments "--timemode 6months --rigidity_start $left --rigidity_end $right --binmerge 2 --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex" Antiproton_LowMakeTimeDependentTemplate 0 0
        cd ..
    done
    cd ..


    mkdir -p 3BartalRotation
    cd 3BartalRotation
    for (( p=0; p<135; p=p+3))
    do
        LeftTimeIndex=$p
        RightTimeIndex=$[p+3]
        mkdir -p $LeftTimeIndex"_"$RightTimeIndex
        cd $LeftTimeIndex"_"$RightTimeIndex
        run_parallel_generic_job --partition c18m --timelimit 24:00 --memory 180 --arguments "--timemode 3BartalRotation --rigidity_start $left --rigidity_end $right --binmerge 2 --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex" Antiproton_LowMakeTimeDependentTemplate 0 0
        cd ..
    done
    cd ..


    mkdir -p 6BartalRotation
    cd 6BartalRotation
    for (( p=0; p<138; p=p+6))
    do
        LeftTimeIndex=$p
        RightTimeIndex=$[p+6]
        mkdir -p $LeftTimeIndex"_"$RightTimeIndex
        cd $LeftTimeIndex"_"$RightTimeIndex
        run_parallel_generic_job --partition c18m --timelimit 24:00 --memory 180 --arguments "--timemode 6BartalRotation --rigidity_start $left --rigidity_end $right --binmerge 2 --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex" Antiproton_LowMakeTimeDependentTemplate 0 0
        cd ..
    done
    cd ..

    cd ..

done

cd ..


