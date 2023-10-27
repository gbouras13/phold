#!/bin/bash

# to make the combined FASTA

input_file="all_deduped_phrogs.fasta"
tsv_file="phrog_index_mapping.tsv"
output_dir="out_phrog_dbs"
database_dir="toy_foldseek_db"
num_threads=8
min_phrog=1
batch_size=10
total_phrog=700

mkdir -p $output_dir

for ((start_phrog = min_phrog; start_phrog <= total_phrog; start_phrog += batch_size)); do
    end_phrog=$((start_phrog + batch_size - 1))
    if ((end_phrog > total_phrog)); then
        end_phrog=$total_phrog
    fi

    output_subdir="${output_dir}/db_${start_phrog}_${end_phrog}"

    echo "Processing phrogs $start_phrog to $end_phrog..."
    
    phold create -i "$input_file" --tsv "$tsv_file" -o "$output_subdir" -t "$num_threads" -d "$database_dir" -f --min_phrog "$start_phrog" --max_phrog "$end_phrog"
done


input_file="all_deduped_phrogs.fasta"
tsv_file="phrog_index_mapping.tsv"
output_dir="out_phrog_dbs"
database_dir="toy_foldseek_db"
num_threads=8
min_phrog=701
batch_size=100
total_phrog=38880

mkdir -p $output_dir

for ((start_phrog = min_phrog; start_phrog <= total_phrog; start_phrog += batch_size)); do
    end_phrog=$((start_phrog + batch_size - 1))
    if ((end_phrog > total_phrog)); then
        end_phrog=$total_phrog
    fi

    output_subdir="${output_dir}/db_${start_phrog}_${end_phrog}"

    echo "Processing phrogs $start_phrog to $end_phrog..."
    
    phold create -i "$input_file" --tsv "$tsv_file" -o "$output_subdir" -t "$num_threads" -d "$database_dir" -f --min_phrog "$start_phrog" --max_phrog "$end_phrog"
done


# to make the foldseek db


cat out_phrog_dbs/*/output3di.fasta > all_phrogs_output3di.fasta
cat out_phrog_dbs/*/outputaa.fasta > all_phrogs_outputaa.fasta


# to make the combined foldseek db

phold createdb --fasta_aa out_phrog_dbs_new/all_phrogs_outputaa.fasta --fasta_3di all_phrogs_output3di.fasta -o PHROG_ProstT5_Foldseek_db_updated -t 16 -f -p all_phrogs





# # to make the foldseek db finally
 
# foldseek concatdbs out_phrog_dbs/db_1_10/foldseek_db/phold_foldseek_database out_phrog_dbs/db_11_20/foldseek_db/phold_foldseek_database  PHROG_ProstT5_Foldseek_db/all_phrog
# foldseek concatdbs out_phrog_dbs/db_1_10/foldseek_db/phold_foldseek_database_ss out_phrog_dbs/db_11_20/foldseek_db/phold_foldseek_database_ss  PHROG_ProstT5_Foldseek_db/all_phrog_ss
# foldseek concatdbs out_phrog_dbs/db_1_10/foldseek_db/phold_foldseek_database_h out_phrog_dbs/db_11_20/foldseek_db/phold_foldseek_database_h  PHROG_ProstT5_Foldseek_db/all_phrog_h
 
# foldseek concatdbs PHROG_ProstT5_Foldseek_db/all_phrog out_phrog_dbs/db_21_30/foldseek_db/phold_foldseek_database  PHROG_ProstT5_Foldseek_db/all_phrog
# foldseek concatdbs PHROG_ProstT5_Foldseek_db/all_phrog_ss out_phrog_dbs/db_21_30/foldseek_db/phold_foldseek_database_ss  PHROG_ProstT5_Foldseek_db/all_phrog_ss
# foldseek concatdbs PHROG_ProstT5_Foldseek_db/all_phrog_h out_phrog_dbs/db_21_30/foldseek_db/phold_foldseek_database_h  PHROG_ProstT5_Foldseek_db/all_phrog_h
 
# # script will loop through the specified range (from 31 to 38880 in increments of 10)
 
# for i in {31..700..10}; do
#     start_range=$i
#     end_range=$((i+9))
 
#     input_db="PHROG_ProstT5_Foldseek_db/all_phrog"
#     output_db="out_phrog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database"
 
#     foldseek concatdbs "$input_db" "$output_db" "$input_db"
 
#     input_db_ss="PHROG_ProstT5_Foldseek_db/all_phrog_ss"
#     output_db_ss="out_phrog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_ss"
 
#     foldseek concatdbs "$input_db_ss" "$output_db_ss" "$input_db_ss"
 
#     input_db_h="PHROG_ProstT5_Foldseek_db/all_phrog_h"
#     output_db_h="out_phrog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_h"
 
#     foldseek concatdbs "$input_db_h" "$output_db_h" "$input_db_h"
# done

# for i in {701..38880..100}; do
#     start_range=$i
#     end_range=$((i+99))
 
#     input_db="PHROG_ProstT5_Foldseek_db/all_phrog"
#     output_db="out_phrog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database"
 
#     foldseek concatdbs "$input_db" "$output_db" "$input_db"
 
#     input_db_ss="PHROG_ProstT5_Foldseek_db/all_phrog_ss"
#     output_db_ss="out_phrog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_ss"
 
#     foldseek concatdbs "$input_db_ss" "$output_db_ss" "$input_db_ss"
 
#     input_db_h="PHROG_ProstT5_Foldseek_db/all_phrog_h"
#     output_db_h="out_phrog_dbs/db_${start_range}_${end_range}/foldseek_db/phold_foldseek_database_h"
 
#     foldseek concatdbs "$input_db_h" "$output_db_h" "$input_db_h"
# done