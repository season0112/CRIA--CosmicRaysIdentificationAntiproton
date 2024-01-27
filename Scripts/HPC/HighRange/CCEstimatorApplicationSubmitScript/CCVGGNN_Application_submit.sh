mkdir -p p0VGGNN/PhyRep2021
cd p0VGGNN/PhyRep2021

bin=("14.1_15.3" "15.3_16.6" "16.6_18" "18_19.5" "19.5_21.1" "21.1_22.8" "22.8_24.7" "24.7_26.7" "26.7_28.8" "28.8_31.1" "31.1_33.5" "33.5_36.1" "36.1_38.9" "38.9_41.9" "41.9_45.1" "45.1_48.5" "48.5_52.2" "52.2_56.1" "56.1_60.3" "60.3_64.8" "64.8_69.7" "69.7_74.9" "74.9_80.5" "80.5_93" "93_108" "108_125" "125_147" "147_175" "175_211" "211_259" "211_250" "250_330" "259_450" "330_525")

for binname in $bin
do
    echo $binname
    mkdir -p $binname
    cd $binname
    run_parallel_generic_job --argument "--rigidityBinOption $binname --issversion PhyRep2021 --ecalBDTCut --trdProtonHeliumCut --physicstriggerCut" --timelimit 24:00 --memory 180 --partition c18m /home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/HighRange/loadmodel_application.py 0 0
    cd ..
done

cd ../..
