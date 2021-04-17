#!/bin/bash
source /reg/g/psdm/etc/psconda.sh
cd /reg/d/psdm/xpp/xpplv9818/results/krapivin/
conda activate ana-4.0.19-py3
mpiexec python3 /reg/d/psdm/xpp/xpplv2818/results/krapivin/process_cube_h5.py $1
