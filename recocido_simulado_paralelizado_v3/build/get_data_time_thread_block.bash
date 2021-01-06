ls -l info_save_test.txt
blockes=(32 64 128 256 )
threads=(1 32 64 85)
for block in "${blockes[@]}"
do
    for thread in "${threads[@]}"
    do
        for repeat in {1..30}
        do 
            line=$(head -n 1 "./save_example_1_1_1/block_$block/save_$thread/$repeat-info_test.txt")
            echo "$line" >> info_save_test_1_1_1.txt
        done
    done
done
echo "---------######## Completado #######----------"
