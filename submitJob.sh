#!/bin/bash
source /reg/g/psdm/etc/psconda.sh
#cd /reg/d/psdm/xpp/xpplv9818/results/krapivin/
cd /cds/data/drpsrcf/xpp/xpplw8919/scratch/krapivin
conda activate ana-4.0.23-py3
#mpiexec python3 /reg/d/psdm/xpp/xpplv2818/results/krapivin/process_cube_h5.py $1
mpiexec python3 /cds/data/drpsrcf/xpp/xpplw8919/scratch/krapivin/process_cube_h5.py --laser_off 1 --pull_from_ffb 1 $1
