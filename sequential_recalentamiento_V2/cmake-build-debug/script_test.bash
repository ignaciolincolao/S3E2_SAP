mkdir "save_example_1_1_1_antigua"
for repeat in {1..30}
do 
    echo "---------######## Estado 15 30 25 ####### En $repeat de 30 ######----------"
    let seed=$(($repeat*$repeat))
    ./sequential 1 1 1 0.97 0.90 10000000000 0.0000000000009 "./save_example_1_1_1_antigua/" "$repeat" $seed
    exit_code=$? 
    while [ $exit_code -gt 0 ];
    do
        ./sequential 1 1 1 0.97 0.90 10000000000 0.0000000000009 "./save_example_1_1_1_antigua/" "$repeat" $seed
        exit_code=$? 
    done
done
echo "---------######## Completado 15 30 25 ####### Completados $repeat de 30 ######----------"
