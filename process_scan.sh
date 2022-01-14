expname=xpplw8419
laserOff=1
ipmLower=10.00
ipm3Lower=-1
detectorThreshold=9.9
relativeScan=1
run=$1

if [ -z "$1" ]
then
    echo "This script is called with 3 parameters, the run number, then the start of the cubing time, then the end of the cubing time."
    exit
fi

mpiexec python process_scan.py --laser_off $laserOff --exp_name $expname --ipm_lower $ipmLower --ipm3_lower $ipm3Lower --detector_threshold $detectorThreshold --relative_scan $relativeScan $run $2
