subdir=$1*
echo "wchen123"
echo "download data from: " $subdir
rsync -av wchen@140.115.34.99:lulin_data/$subdir ./
