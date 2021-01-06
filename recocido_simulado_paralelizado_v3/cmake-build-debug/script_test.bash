blockes=(32 64 128 256 )
threads=(1 32 64 85)
mkdir "save_example_15_30_25"
for repeat in {1..30}
do 
    for block in "${blockes[@]}"
    do
        mkdir "./save_example_15_30_25/block_$block"
        for thread in "${threads[@]}"
        do
            mkdir "./save_example_15_30_25/block_$block/save_$thread"
            echo "---------######## Estado 15 30 25 ####### En $repeat de 30 ######----------"
            let seed=$(($repeat+$block*$thread))
            ./paralelizado 15 30 25 0.98 0.90 100000 0.00000009 $block $thread "./save_example_15_30_25/block_$block/save_$thread/" "$repeat" $seed
            exit_code=$? 
            while [ $exit_code -gt 0 ];
            do
                ./paralelizado 15 30 25 0.98 0.90 100000 0.00000009 $block $thread "./save_example_15_30_25/block_$block/save_$thread/" "$repeat" $seed
                exit_code=$? 
            done

        done
    done
done
echo "---------######## Completado 15 30 25 ####### Completados $repeat de 5 ######----------"
echo "---------### bloques 16 32 64 128 256  #### hilos 1 16 32 64 85 #####----------" 

