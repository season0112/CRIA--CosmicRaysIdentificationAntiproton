mkdir -p Submit_TemplateFitForSyematicEff_Dependent_Low
cd Submit_TemplateFitForSyematicEff_Dependent_Low

    
mkdir -p 6months
cd 6months
for (( p=0; p<21; p=p+1))
do
    LeftTimeIndex=$p
    RightTimeIndex=$[p+1] 
    mkdir -p $LeftTimeIndex"_"$RightTimeIndex
    cd $LeftTimeIndex"_"$RightTimeIndex
    run_parallel_generic_job --partition c18m --timelimit 24:00 --memory 180 --arguments "--timemode 6months --rigidity_start 2 --rigidity_end 17 --binmerge 2 --ParametrilizedMode No --GenerateNumbers 1 --SysErrEffStudyMode Yes --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex" Antiproton_LowTF1D_TimeStamp_3and6_BartalRotations 0 0
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
    run_parallel_generic_job --partition c18m --timelimit 24:00 --memory 180 --arguments "--timemode 3BartalRotation --rigidity_start 2 --rigidity_end 17 --binmerge 2 --ParametrilizedMode No --GenerateNumbers 1 --SysErrEffStudyMode Yes --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex" Antiproton_LowTF1D_TimeStamp_3and6_BartalRotations 0 0
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
    run_parallel_generic_job --partition c18m --timelimit 24:00 --memory 180 --arguments "--timemode 6BartalRotation --rigidity_start 2 --rigidity_end 17 --binmerge 2 --ParametrilizedMode No --GenerateNumbers 1 --SysErrEffStudyMode Yes --StartTimeIndex $LeftTimeIndex --EndTimeIndex $RightTimeIndex" Antiproton_LowTF1D_TimeStamp_3and6_BartalRotations 0 0
    cd ..
done
cd ..


cd ..


