mkdir -p Submit_Antiproton_Intermediate_MakeTimeAveragedTemplate
cd Submit_Antiproton_Intermediate_MakeTimeAveragedTemplate

#binmerge_all=("1")
#binmerge_all=("2")
binmerge_all=("1" "2")

for binmerge in $binmerge_all
do
    mkdir -p "binmerge"$binmerge
    cd "binmerge"$binmerge

    for (( i=10; i<30; i=i+binmerge))
    do
        echo Now is: rigidity_start = $i ,  rigidity_end = $[i+binmerge]

        left=$i
        right=$[i+binmerge]

        mkdir -p $left"_"$right
        cd $left"_"$right

        mkdir -p 2016paper
        cd 2016paper
        run_parallel_generic_job --timelimit 24:00 --memory 180 --arguments "--richcut free --trackerpattern 0124 --issversion 2016paper --rigidity_start $left --rigidity_end $right --binmerge $binmerge --signalefficiency 0.8 --Testmode Yes --FullRange Yes" Antiproton_IntermediateMakeTemplate 0 0
        cd ..

        mkdir -p PhysicsReport
        cd PhysicsReport
        run_parallel_generic_job --timelimit 24:00 --memory 180 --arguments "--richcut free --trackerpattern 0124 --issversion PhysicsReport --rigidity_start $left --rigidity_end $right --binmerge $binmerge --signalefficiency 0.8 --Testmode Yes --FullRange Yes" Antiproton_IntermediateMakeTemplate 0 0
        cd ..

        mkdir -p pass7.8
        cd pass7.8
        run_parallel_generic_job --timelimit 24:00 --memory 180 --arguments "--richcut free --trackerpattern 0124 --issversion pass7.8 --rigidity_start $left --rigidity_end $right --binmerge $binmerge --signalefficiency 0.8 --Testmode Yes --FullRange Yes" Antiproton_IntermediateMakeTemplate 0 0
        cd ..

        cd ..

    done

    cd ..

done

cd ..


