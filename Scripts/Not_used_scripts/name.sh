for file in ExampleAnalysis_*
do
    name=${file}
    name2=$(echo $name | sed 's/00000/00015/g')
    echo ${name2}
    mv ${file} ${name2}
done
