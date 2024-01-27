mkdir -p Submit_Antiproton_Intermediate_MakeTimeDependentTemplate
cd Submit_Antiproton_Intermediate_MakeTimeDependentTemplate

for (( i=10; i<30; i=i+2))
do
    echo Now is: rigidity_start = $i ,  rigidity_end = $[i+1]

    left=$i
    right=$[i+1]

    mkdir -p $left"_"$right
    cd $left"_"$right
   
    ''' 
    mkdir -p 6months
    cd 6months
    for (( p=0; p<21; p=p+1))
    do
        LeftTimeIndex=$p
        RightTimeIndex=$[p+1]
        mkdir -p $LeftTimeIndex"_"$RightTimeIndex
        cd $LeftTimeIndex"_"$RightTimeIndex
        run_parallel_generic_job  --timelimit 24:00 --memory 180 --arguments "--richcut free --trackerpattern 0124 --rigidity_start $left --rigidity_end $right --binmerge 2 --signalefficiency 0.8 --timemode 6months --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex --FullRange Yes" Antiproton_IntermediateMakeTemplate_TimeDependent 0 0
        cd ..
    done
    cd ..
    '''

    '''
    mkdir -p 3BartalRotation
    cd 3BartalRotation
    for (( p=0; p<135; p=p+3))
    do
        LeftTimeIndex=$p
        RightTimeIndex=$[p+3]
        mkdir -p $LeftTimeIndex"_"$RightTimeIndex
        cd $LeftTimeIndex"_"$RightTimeIndex
        run_parallel_generic_job  --timelimit 24:00 --memory 180 --arguments "--richcut free --trackerpattern 0124 --rigidity_start $left --rigidity_end $right --binmerge 2 --signalefficiency 0.8 --timemode 3BartalRotation --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex --FullRange Yes" Antiproton_IntermediateMakeTemplate_TimeDependent 0 0
        cd ..
    done
    cd ..
    '''

    mkdir -p 6BartalRotation
    cd 6BartalRotation
    for (( p=0; p<138; p=p+6))
    do
        LeftTimeIndex=$p
        RightTimeIndex=$[p+6]
        mkdir -p $LeftTimeIndex"_"$RightTimeIndex
        cd $LeftTimeIndex"_"$RightTimeIndex
        run_parallel_generic_job  --timelimit 24:00 --memory 180 --arguments "--richcut free --trackerpattern 0124 --rigidity_start $left --rigidity_end $right --binmerge 2 --signalefficiency 0.8 --timemode 6BartalRotation --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex --FullRange Yes" Antiproton_IntermediateMakeTemplate_TimeDependent 0 0
        cd ..
    done
    cd ..

    cd ..

done

cd ..


