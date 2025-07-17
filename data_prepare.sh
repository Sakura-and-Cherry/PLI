file_path="/mnt/Database/Datasets/SOSD/"
output_path="/mnt/Database/Datasets/SOSD/pgm/"

# "books_800M_uint64" "osm_cellids_800M_uint64"
# Build the PGM index for each dataset
for dataset in "books_200M_uint32"
do
  dataset_name=$file_path$dataset
  unique_dataset_name=$file_path$dataset"_unique"
  echo "Preparing dataset: $dataset_name"
  ./build/data_prepare_uint32 $dataset_name f
  ./build/data_prepare_uint32 $unique_dataset_name t
done

for dataset in "fb_200M_uint64" "wiki_ts_200M_uint64" "books_800M_uint64" "osm_cellids_800M_uint64"
do
  dataset_name=$file_path$dataset
  unique_dataset_name=$file_path$dataset"_unique"
  echo "Preparing dataset: $dataset_name"
  ./build/data_prepare_uint64 $dataset_name f
  ./build/data_prepare_uint64 $unique_dataset_name t
done