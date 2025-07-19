file_path="/mnt/Database/Datasets/SOSD/"
output_path="/mnt/Database/Datasets/SOSD/pgm/"

mkdir -p $output_path
mkdir -p "./log/"
# Build the PGM index for each dataset

partition_type="optimal"
#"0.0001" "0.0005" "0.001" "0.005"
for dataset in "test"
  do
    log_path="./log/$dataset""_unique_""$partition_type/build.log"
    mkdir -p "./log/$dataset""_unique_""$partition_type"
    rm -f $log_path
    for lambda in "0.0008" "0.00085" "0.0009"
    do
      dataset_name=$file_path$dataset"_unique"
      index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
      echo "Processing dataset: $dataset_name"
      ./build/pgm_build_uint32 $dataset_name $index_name f $lambda $partition_type $log_path
      echo " "
  done
  echo "——————————————"
done
