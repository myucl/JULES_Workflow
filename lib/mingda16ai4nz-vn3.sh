#!/bin/bash

# bash files are used on LINUX systems.
#
# Author: Carolina Duran Rojas
#
# Created: Jun 2023
#
# Last modification date: 10 July 2021
#
# This script is used to:
# . access the suite to run JULES 
# . run JULES
# . stop JULES runs  
#
# Usage: ./ai4nz.sh

#set -x

# modules to run suites



module load oneapi/compilers/24.2.0
module load oneapi/mpi/24.2.0
module load netcdf/intel2024.2.0/4.9.2
module load netcdf/intel2024.2.0/fortran/4.6.1

# libraries to run suites

#export HDF5_LIBDIR=/apps/jasmin/supported/libs/hdf5/intel2024.2.0/1.14.4-2/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_LIBDIR
#export NETCDF_FORTRAN_ROOT=/apps/jasmin/supported/libs/netcdf/intel2024.2.0/fortran/4.6.1
export HDF5_LIBDIR=/apps/jasmin/supported/libs/hdf5/intel2024.2.0/1.14.4-2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_LIBDIR
export NETCDF_FORTRAN_ROOT=/apps/jasmin/supported/libs/netcdf/intel2024.2.0/fortran/4.6.1
# remove files from prevoius runs
#OUTPUT_FOLDER=/gws/nopw/j04/uknetzero/mingda16/u-cx502-species_calibration/ENSEMBLE/
#uncomment when all ready
#$( rm -rf /work/scratch-pw3/mingda16/u-cx502-species_calibration/ENSEMBLE/* )


OUTPUT_FOLDER=/work/scratch-pw3/mingda16/u-cx502-species_calibration/ENSEMBLE/

#$( rm -rf /gws/nopw/j04/uknetzero/public/jules/ENSEMBLE/* )


# number of points to model
ENSEMBLE_N=`wc -l  < coordinates.csv`
ENSEMBLE_NUM=$(( $ENSEMBLE_N -1 ))

# username
USERNAME=`whoami`

echo `date`

# path of the suite
cd /home/users/mingda16/roses/species_calibration

# changing the  rose-suite.conf file for the suite
sed -i 's/ENSEMBLE_NUM=[0-9]*/ENSEMBLE_NUM='$ENSEMBLE_NUM'/g' rose-suite.conf

# run suite with the GUI
#rose suite-run --new 
# run suite without the GUI

cylc vip species_calibration

echo `date`

#status_file=$HOME/cylc-run/species_calibration/run1/log/job/19610101T0000+0100/JULES_000001/NN/job.status
status_file=$HOME/cylc-run/species_calibration/run1/log/job/*/JULES_000001/NN/job.status

# Wait until the status file is generated
while [[ ! -f $status_file ]]; do
  sleep 1m
done

for (( ens_num = 1; ens_num <= $ENSEMBLE_NUM; ens_num++ )); do
  while : ; do
    run=$( cylc ping -v --host=cylc1.jasmin.ac.uk species_calibration > running.out )
    running=$( grep -c "Running" running.out )
    if [[ $running -eq 1 ]]; then 
      sleep 2m
      echo `date`
    else
      stalled=$(cylc log -f a species_calibration | grep -c stalled)
      if [[ $stalled -gt 0 ]]; then
        exit
      else
         break
      fi
    fi
  done < running.out 
done

echo `date`

cylc stop --now --now species_calibration/run* ; cylc clean species_calibration

#cp /work/scratch-pw3/mingda16/u-cx502-species_calibration/ENSEMBLE/#RAll_UF_N000100_S000000001/000*/JULES7.0-RED1.1_*_NT.monthly.nc  $HOME/jules_output



echo `date`
exit
 
