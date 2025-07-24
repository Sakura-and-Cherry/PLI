#!/bin/bash
file_path="data/SOSD/"
output_path="data/SOSD/pgm/"


mkdir -p $output_path

partition_type="uniform"

for dataset in "books_200M_uint32"
  do
    log_path="./log/$dataset""_unique_""$partition_type/query.log"
    mkdir -p "./log/$dataset""_unique_""$partition_type"
    rm -f $log_path
    for lambda in "0.00075" "0.0008"
      do
      dataset_name=$file_path$dataset"_unique"
      index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
      echo "Querying index: $index_name lambda: $lambda"
      ./build/pgm_query_uint32 $dataset_name $index_name $lambda $partition_type $log_path
      echo " "
    done
  echo "——————————————"
done