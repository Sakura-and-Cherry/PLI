file_path="/mnt/Database/Datasets/SOSD/"
output_path="/mnt/Database/Datasets/SOSD/pgm/"

mkdir -p $output_path
# Build the PGM index for each dataset

partition_type="greedy"

#"fb_200M_uint64" "wiki_ts_200M_uint64" "books_800M_uint64" "osm_cellids_800M_uint64"
for dataset in "osm_cellids_800M_uint64"
do
   for lambda in "0.0001" "0.0005" "0.001" "0.005" "0.05" "0.5"
  do
    dataset_name=$file_path$dataset"_unique"
    index_name=$output_path$dataset"_unique_"$partition_type\_$lambda".idx"
    log_path="./log/$dataset""_unique_""$partition_type/build.log"
    mkdir -p "./log/$dataset""_unique_""$partition_type"
    touch $log_path
    echo "Processing dataset: $dataset_name lambda: $lambda"
    ./build/pgm_build_uint64 $dataset_name $index_name f $lambda $partition_type $log_path
    echo " "
  done
  echo "——————————————"
done