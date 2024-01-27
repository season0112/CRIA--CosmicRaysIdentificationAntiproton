
cd /hpcwork/jara0052/sichen/Measuringtime/TimeDependent

Degree="25"
SafetyFactor="1.2"
UnixtimeEdge_6B=("1305417600" "1319414400" "1333411200" "1347408000" "1361404800" "1375401600" "1389398400" "1403395200" "1417392000" "1431388800" "1445385600" "1459382400" "1473379200" "1487376000" "1501372800" "1515369600" "1529366400" "1543363200" "1557360000" "1571356800" "1585353600" "1599350400" "1613347200" "1627344000")

for i in {1..23}
do
    cd ${UnixtimeEdge_6B[$i]}_${UnixtimeEdge_6B[$i+1]}
    ac_merge --input /hpcwork/jara0052/sichen/Measuringtime/TimeDependent/${UnixtimeEdge_6B[$i]}_${UnixtimeEdge_6B[$i+1]}/MeasuringTime*.root --resultprefix MeasuringTime_pass7.8_11_2021_GEOMETRIC${Degree}_${SafetyFactor}_${UnixtimeEdge_6B[$i]}_${UnixtimeEdge_6B[$i+1]} --resultdir /hpcwork/jara0052/sichen/Measuringtime/TimeDependent/
    cd ..
done
cd /home/bo791269/Software/AntiprotonAnalysis/Scripts/SubmitScript


