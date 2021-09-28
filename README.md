# PSANA Binning Code for the XPP Hutch
Here I have written some binning code for LCLS experiments in the XPP hutch.

This code here can process the full detector image and output a code for a particular run. Generally, this process is significantly faster than the facility provided cube code as it directly interfaces with psana.

## How to run this code
You need to first the psana anaconda environment to your path. You can do this by running (or adding this line at the end of your .bashrc)
```
source /reg/g/psdm/etc/psconda.sh
```
Then to run the code
```
ssh psana
cd <directory where code is located>
conda activate ana-4.0.22-py3 # this sets your python to python 3 instead of the default python 2. Here I have provided the version of the environment I was using, but generally this code should work for the current version
mpiexec python3 process_phiscan_h5.py #runnumberhere
```
You may also modify submitJob_phiscan.sh appropriately (changing directories to those that correspond to your directory). and then run
```
sbatch -n 20 --exclusive --ntasks-per-node=1 submitJob_phiscan_lasoff.sh #runnumberhere
```
The above code is best to use if you need critically fast processing time during an experiment. In post beamtime processing, I suggest using only 6 nodes without the exclusive or the ntasks-per-node option above. These last two options force the process to be able to use the entire memory contents of a node.

If you have questions, you can reach out to krapivin at [inset commonly used name of Leland Stanford Jr. University].edu

## Future plans
1. More command line options
2. Additional batch scripts to quickly run the cube code with certain parameters set
3. Removal of all hard coded parameters in the code to make it more generic and remove requirements to edit the code during a beamtime.