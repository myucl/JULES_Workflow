#!/bin/bash

# bash files are used on LINUX systems.
#
# Author: Carolina Duran Rojas
#
# Created: Jun 2023
#
# Last modification date: 7 August 2021
#
# This script is used to 
# . send a file of coordinates to JASMIN 
# . connect to JASMIN
# . run JULES on JASMIN
# . retrieve nc files
#
# Usage: ./jas-jls.sh coordinate_file 
#
# Optional Usage: ./jas-jls.sh coordinate_file > /dev/null

#set -x

here=`pwd`
#echo $here
#echo `date`

# uncomment this line if you give csv file from the command line
# reading coordinates file from the command line
coord_file=$1

# this step is user dependent, due to credentials
# ?? Uncomment following line, when someone else uses this script. 
#val $(ssh-agent -s); ssh-add ~/.ssh/id_rsa_jasmin

# retrieving JASMIN username 
# ?? Uncomment the following 2 lines 
#echo "JASMIN username: "
#read username
# ?? Comment this when prevoius lines are uncommented
USERNAME=mingda16

# transfer original coordinates file to JASMIN
# ?? Comment the following 3 lines. 
#echo "coordinates file to JASMIN"
#read coord_file
scp $coord_file $USERNAME@xfer-vm-01.jasmin.ac.uk:~/coordinates.csv
#rsync -ravzpu $coord_file $USERNAME@xfer-vm-01.jasmin.ac.uk:~/coordinates.csv


# ?? Comment following line when previous three are uncommented
#rsync -arpuvz coord.csv caroduro@xfer3.jasmin.ac.uk:/home/users/$username/coordinates.csv

# count the number of lines for the coordinates file
ENSEMBLE_N=`wc -l  < $coord_file`
#ENSEMBLE_N=`wc -l  < coord.csv`
ENSEMBLE_NUM=$(( $ENSEMBLE_N-1 ))

echo `date`
# connect to the JASMIN HPC and run the bash file (run the suite to have output files) 
#ssh -t -AX $USERNAME@cylc.jasmin.ac.uk " ./ai4nz.sh; bash -l "
#ssh -t -AX $USERNAME@cylc.jasmin.ac.uk "nohup  ./ai4nz-vn2.sh"
#ssh -t -AX $USERNAME@cylc.jasmin.ac.uk "./mingda16ai4nz-vn2.sh >/dev/null "
#ssh -t -AX $USERNAME@cylc.jasmin.ac.uk "./mingda16ai4nz-vn3.sh"
#ssh -t username@jasmin "bash --login myscript.sh"
ssh -t -AX $USERNAME@cylc2 "bash --login ./mingda16ai4nz-vn3.sh  " 

#ps_status=`ps -u`
#OUTPUT_FOLDER='/gws/nopw/j04/uknetzero/$USERNAME/u-cx502-species_calibration/ENSEMBLE/'
OUTPUT_FOLDER='/work/scratch-pw3/mingda16/u-cx502-species_calibration/ENSEMBLE/'
# transfering back JULES outputs
#for (( ens_num = 1; ens_num <= $ENSEMBLE_NUM; ens_num++ )); do
  #scp $USERNAME@xfer3.jasmin.ac.uk:$OUTPUT_FOLDER/JULES7.0-RED1.1_rcp26_06_*.monthly.nc $here

#OUTPUT_FOLDER='/home/users/mingda16/jules_output/'
# transfering back JULES outputs
#for (( ens_num = 1; ens_num <= $ENSEMBLE_NUM; ens_num++ )); do
  #scp $USERNAME@xfer3.jasmin.ac.uk:$OUTPUT_FOLDER/JULES7.0-RED1.1*.monthly.nc $here
#done
echo `date`

scp -r mingda16@cylc2:/work/scratch-pw3/mingda16/u-cx502-species_calibration/ENSEMBLE/*/*/JULES7.0-RED1.1_rcp85_15_JULESAPPNEW*NT.monthly.nc  $here/
#change the file name without ALLYCT1 when ready