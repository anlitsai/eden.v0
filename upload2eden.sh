folder=$1

python edenTransfer/store_fits.py --dryrun /home/altsai/eden/$folder/ |tee $folder'_dryrun.log'
python edenTransfer/store_fits.py --verbose /home/altsai/eden/$folder/ |tee $folder'_2eden.log'

