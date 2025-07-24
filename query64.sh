#!/bin/bash
file_path="/home/andy/Projects/Dataset/SOSD"
output_path="/home/andy/Projects/Dataset/SOSD/pgm"



mkdir -p $output_path

partition_type="greedy"
# for dataset in "fb_200M_uint64"
# do
#   log_path="./log/$dataset""_unique_""$partition_type/query.log"
#   mkdir -p "./log/$dataset""_unique_""$partition_type"
#   rm -f $log_path
#    for lambda in "0.0001" "0.0002" "0.0003" "0.00035" "0.000375" "0.0004" "0.0005"
#   do
#     dataset_name=$file_path$dataset"_unique"
#     index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
#     echo "Querying index: $index_name lambda: $lambda"
#     ./build/pgm_query_uint64 $dataset_name $index_name $lambda $partition_type $log_path
#     echo " "
#   done
#   echo "——————————————"
# done


for dataset in "books_800M_uint64"
do
  log_path="./log/$dataset""_unique_""$partition_type/query.log"
  mkdir -p "./log/$dataset""_unique_""$partition_type"
  rm -f $log_path
   for lambda in "0.0004" "0.0005" "0.00055" "0.0006" "0.000625" "0.00065" "0.0007"
  do
    dataset_name=$file_path$dataset"_unique"
    index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
    echo "Querying index: $index_name lambda: $lambda"
    ./build/pgm_query_uint64 $dataset_name $index_name $lambda $partition_type $log_path
    echo " "
  done
  echo "——————————————"
done