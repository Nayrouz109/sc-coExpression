


for file in /space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/CMC/Raw/*.rds; do
    base_name=$(basename $file)
    new_name=$(echo $base_name | sed 's/CMC_SZ/MSSM-2024/')
    cp "$file" /space/scratch/nairuz-rewiring/data/0.processed-raw/$new_name
done



for file in /space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/SZBDMulti-Seq/Raw/*.rds; do
    base_name=$(basename $file)
    new_name=$(echo $base_name | sed 's/SZBDMulti-Seq_SZ/SZBDMulti-2024/')
    cp "$file" /space/scratch/nairuz-rewiring/data/0.processed-raw/$new_name
done




for file in /space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/CMC/Raw/*.rds; do
    base_name=$(basename $file)
    new_name=$(echo $base_name | sed 's/CMC_SZ/MSSM-2024/')
    cp "$file" /cosmos/data/downloaded-data/AD-meta/ruzika-2024/data/$new_name
done



for file in /space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/SZBDMulti-Seq/Raw/*.rds; do
    base_name=$(basename $file)
    new_name=$(echo $base_name | sed 's/SZBDMulti-Seq_SZ/SZBDMulti-2024/')
    cp "$file" /cosmos/data/downloaded-data/AD-meta/ruzika-2024/data/$new_name
done


