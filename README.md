# psana_binning_code
some rough binning code LCLS psana XPP experiments
will attempt to make this more user friendly for later experiments.

## getting it to run
You need to first the psana anaconda environment to your path. You can do this by running (or adding this line at the end of your .bashrc)
```
source /reg/g/psdm/etc/psconda.sh
```
Then to run the code
```
ssh psana
cd <directory where code is located>
conda activate ana-4.0.22-py3 # this sets your python to python 3 instead of the default python 2.
mpiexec python3 process_phiscan_h5.py #runnumberhere
```
You may also modify submitJob_phiscan.sh appropriately (changing directories to those that correspond to your directory). and then run
```
sbatch -n 20 --exclusive --ntasks-per-node=1 submitJob_phiscan_lasoff.sh #runnumberhere
```

If you have questions, you can reach out to krapivin at [inset commonly used name of Leland Stanford Jr. University].edu
