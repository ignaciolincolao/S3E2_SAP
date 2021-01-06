mkdir "save_example_15_30_25_v2"
for repeat in {1..30}
do 
    echo "---------######## Estado 1 1 1 ####### En $repeat de 30 ######----------"
    let seed=$(($repeat*$repeat))
    ./sequential 15 30 25 0.97 0.90 1 0.0000000000009 "./save_example_15_30_25_v2/" "$repeat" $seed
    exit_code=$? 
    while [ $exit_code -gt 0 ];
    do
        ./sequential 15 30 25 0.90 0.90 1 0.0000000000000009 "./save_example_15_30_25_v2/" "$repeat" $seed
        exit_code=$? 
    done
done
echo "---------######## Completado 1 1 1 ####### Completados $repeat de 30 ######----------"
