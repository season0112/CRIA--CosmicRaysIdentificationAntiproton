patterrnlist=("-1" "0" "1" "2" "3" "4" "5")
#issversion=("published2016" "PhyRep2021" "pass7.8")
issversion=("published2016" "PhyRep2021")

for issused in $issversion
do
    for pattern in $patterrnlist
    do
        Unfolding_Bayes --binningversion 525version --issversion $issused --pattern $pattern --ifVGGNN No
    done
    Unfolding_Bayes --binningversion 525version --issversion $issused --pattern 0 --ifVGGNN Yes
done

