#!/bin/bash
file_path="/home/andy/Projects/Dataset/SOSD/"
output_path="/home/andy/Projects/Dataset/SOSD/pgm/"



mkdir -p $output_path
# Build the PGM index for each dataset



#"fb_200M_uint64"
# for dataset in "fb_200M_uint64"
# do
#   log_path="./log/$dataset""_unique_""$partition_type/build.log"
#   mkdir -p "./log/$dataset""_unique_""$partition_type"
#   rm -f $log_path
#    for lambda in "0.0001" "0.0002" "0.0003" "0.00035" "0.000375" "0.0004" "0.0005"
#    do
#     dataset_name=$file_path$dataset"_unique"
#     index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
#     echo "Processing dataset: $dataset_name lambda: $lambda"
#     ./build/pgm_build_uint64 $dataset_name $index_name f $lambda $partition_type $log_path
#     echo " "
#   done
#   echo "——————————————"
# done

#"books_800M_uint64" greedy
partition_type="optimal"
for dataset in "fb_10M_uint64"
do
  log_path="./log/$dataset""_unique_""$partition_type/build.log"
  mkdir -p "./log/$dataset""_unique_""$partition_type"
  rm -f $log_path
   for lambda in "0.01" "0.02" "0.03" "0.04" "0.05" "0.06" "0.07"
   do
    dataset_name=$file_path$dataset"_unique"
    index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
    echo "Processing dataset: $dataset_name lambda: $lambda"
    ./build/pgm_build_uint64 $dataset_name $index_name f $lambda $partition_type $log_path
    echo " "
  done
  echo "——————————————"
done


# "osm_cellids_800M_uint64"
# for dataset in "osm_cellids_800M_uint64"
# do
#   log_path="./log/$dataset""_unique_""$partition_type/build.log"
#   mkdir -p "./log/$dataset""_unique_""$partition_type"
#   rm -f $log_path
#    for lambda in "0.0001" "0.0002" "0.0003" "0.00035" "0.000375" "0.0004" "0.0005"
#    do
#     dataset_name=$file_path$dataset"_unique"
#     index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
#     echo "Processing dataset: $dataset_name lambda: $lambda"
#     ./build/pgm_build_uint64 $dataset_name $index_name f $lambda $partition_type $log_path
#     echo " "
#   done
#   echo "——————————————"
# done
