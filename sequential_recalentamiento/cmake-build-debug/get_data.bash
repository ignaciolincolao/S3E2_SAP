alphas="15_30_25"
ls -l info_save_test_$alphas.txt
for repeat in {1..30}
    do 
        line=$(head -n 1 "./save_example_15_30_25_v3/$repeat-info_test.txt")
        echo "$line" >> info_save_test_v3_$alphas.txt
    done
echo "---------######## Completado #######----------"
