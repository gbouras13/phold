#!/bin/bash

# to make the combined FASTA

# file from megapharokka
# no need for tsv
input_file="phold_envhog_db/EnVhog_consensus_renamed.fasta"
output_dir="out_envhog_dbs"
database_dir="toy_foldseek_db"
num_threads=8
min_envhog=1
batch_size=1000
total_envhog=623000

#commit 1472f4f

mkdir -p $output_dir

for ((start_envhog = min_envhog; start_envhog <= total_envhog; start_envhog += batch_size)); do
    end_envhog=$((start_envhog + batch_size - 1))
    if ((end_envhog > total_envhog)); then
        end_envhog=$total_envhog
    fi

    output_subdir="${output_dir}/db_${start_envhog}_${end_envhog}"

    echo "Processing envhogs $start_envhog to $end_envhog..."
    
    phold create -i "$input_file" --tsv "$tsv_file" -o "$output_subdir" -t "$num_threads" -d "$database_dir" -f --envhog_flag --envhog_start "$start_envhog" --envhog_batch_size "$batch_size"
done



# to make the foldseek db finally
 
# foldseek concatdbs out_envhog_dbs/db_1_10/foldseek_db/phold_foldseek_database out_envhog_dbs/db_11_20/foldseek_db/phold_foldseek_database  envhog_ProstT5_Foldseek_envhog_db/all_envhog
# foldseek concatdbs out_envhog_dbs/db_1_10/foldseek_db/phold_foldseek_database_ss out_envhog_dbs/db_11_20/foldseek_db/phold_foldseek_database_ss  envhog_ProstT5_Foldseek_envhog_db/all_envhog_ss
# foldseek concatdbs out_envhog_dbs/db_1_10/foldseek_db/phold_foldseek_database_h out_envhog_dbs/db_11_20/foldseek_db/phold_foldseek_database_h  envhog_ProstT5_Foldseek_envhog_db/all_envhog_h
 
# foldseek concatdbs envhog_ProstT5_Foldseek_envhog_db/all_envhog out_envhog_dbs/db_21_30/foldseek_db/phold_foldseek_database  envhog_ProstT5_Foldseek_envhog_db/all_envhog
# foldseek concatdbs envhog_ProstT5_Foldseek_envhog_db/all_envhog_ss out_envhog_dbs/db_21_30/foldseek_db/phold_foldseek_database_ss  envhog_ProstT5_Foldseek_envhog_db/all_envhog_ss
# foldseek concatdbs envhog_ProstT5_Foldseek_envhog_db/all_envhog_h out_envhog_dbs/db_21_30/foldseek_db/phold_foldseek_database_h  envhog_ProstT5_Foldseek_envhog_db/all_envhog_h
 
# # script will loop through the specified range (from 31 to 38880 in increments of 10)
 
# for i in {31..700..10}; do
#     start_range=$i
#     end_range=$((i+9))
 
#     input_db="envhog_ProstT5_Foldseek_envhog_db/all_envhog"
#     output_db="out_envhog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database"
 
#     foldseek concatdbs "$input_db" "$output_db" "$input_db"
 
#     input_db_ss="envhog_ProstT5_Foldseek_envhog_db/all_envhog_ss"
#     output_db_ss="out_envhog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_ss"
 
#     foldseek concatdbs "$input_db_ss" "$output_db_ss" "$input_db_ss"
 
#     input_db_h="envhog_ProstT5_Foldseek_envhog_db/all_envhog_h"
#     output_db_h="out_envhog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_h"
 
#     foldseek concatdbs "$input_db_h" "$output_db_h" "$input_db_h"
# done

# for i in {701..38880..100}; do
#     start_range=$i
#     end_range=$((i+99))
 
#     input_db="envhog_ProstT5_Foldseek_envhog_db/all_envhog"
#     output_db="out_envhog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database"
 
#     foldseek concatdbs "$input_db" "$output_db" "$input_db"
 
#     input_db_ss="envhog_ProstT5_Foldseek_envhog_db/all_envhog_ss"
#     output_db_ss="out_envhog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_ss"
 
#     foldseek concatdbs "$input_db_ss" "$output_db_ss" "$input_db_ss"
 
#     input_db_h="envhog_ProstT5_Foldseek_envhog_db/all_envhog_h"
#     output_db_h="out_envhog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_h"
 
#     foldseek concatdbs "$input_db_h" "$output_db_h" "$input_db_h"
# done