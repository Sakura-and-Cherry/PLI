
file_path="/mnt/Database/Datasets/SOSD/"
output_path="/mnt/Database/Datasets/SOSD/pgm/"

mkdir -p $output_path

partition_type="greedy"

for dataset in "books_200M_uint32"
  do
    for lambda in "0.0005" "0.0006" "0.0007" "0.0008" "0.0009" "0.001"
      do
      dataset_name=$file_path$dataset"_unique"
      index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
      log_path="./log/$dataset""_unique_""$partition_type/query.log"
      mkdir -p "./log/$dataset""_unique_""$partition_type"
      touch $log_path
      echo "Querying index: $index_name lambda: $lambda"
      ./build/pgm_query_uint32 $dataset_name $index_name $lambda $partition_type $log_path
      echo " "
    done
  echo "——————————————"
done