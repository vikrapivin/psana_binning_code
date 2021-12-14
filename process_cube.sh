expname=xpplw8419
laserOff=1
ipmLower=10.00
ttamp=0.005
numberOfBins=400
detectorThreshold=10
timeStart=$2
timeEnd=$3
run=$1
ttCalibration=0.00156375
ttOffset=0.88457732

if [ -z "$1" ]
then
    echo "This script is called with 3 parameters, the run number, then the start of the cubing time, then the end of the cubing time."
    exit
fi

mpiexec python process_cube_h5.py --laser_off $laserOff --bin_number $numberOfBins --time_low $timeStart --time_high $timeEnd --exp_name $expname --ipm_lower $ipmLower --tt_amp $ttamp --detector_threshold $detectorThreshold --tt_offset $ttOffset --tt_calibration $ttCalibration $run $4
