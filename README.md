# PSANA Binning Code for the XPP Hutch
Here I have written some binning code for LCLS experiments in the XPP hutch.

This code here can process the full detector image and output a code for a particular run. Generally, this process is significantly faster than the facility provided cube code as it directly interfaces with psana.

## How to run this code
Please first log in with your SLAC account to pslogin.
For the first time you run this code, you need to add the psana anaconda environment to your path. You can do this by running (or adding this line at the end of your .bashrc, credits Vincent)
```
#sit_setup ana-current
if [ $(echo $HOSTNAME | grep -ic -e "psana" -e "mon0" -e "drp-srcf") -eq 1 ] 
then
        #echo "on a psana node"
        source /reg/g/psdm/etc/psconda.sh -py3 &> /dev/null
elif [ $(echo $HOSTNAME | grep -ic -e "ioc" -e "control" -e "monitor") -eq 1 ]
then
        echo "on an ioc or non DAQ control room node."
        source /reg/g/pcds/setup/epicsenv-cur.sh
        export PSPKG_ROOT=/reg/g/pcds/pkg_mgr
        source /reg/g/pcds/setup/pcds_shortcuts.sh
elif [ $(echo $HOSTNAME | grep -ic -e "psbuild" -e "psdev") -eq 1 ]
then
        source pcds_conda
        source /reg/g/pcds/setup/plcenv.sh
        source /reg/g/pcds/setup/epicsenv-cur.sh
        source /reg/g/pcds/setup/pcds_shortcuts.sh
else
        export PATH="$PATH:/reg/common/tools/bin"
fi
```
Then to run the code
```
ssh psana
cd <directory where code is located>
mpiexec python3 process_cube.sh #runnumberhere #timeStart #timeEnd #"--num_events maxNumberOfEventsHere"
```
Note: you may need to run something like conda activate `ana-4.0.35-py3`. This sets your python to python 3 instead of the default python 2. Here I have provided the version of the environment I was using, but generally this code should work for the current version.

The above code runs on the psanagpu node you are currently located on. Typically this is not how you want to run cubing code. This way of running the code is useful to test this code. You should instead run something like:

```
sbatch -n 20 --exclusive -p psanaq process_cube.sh 355 -5 20
```
to do a time cube and run
```
sbatch -n 4 --exclusive -p psanaq process_scan.sh 23
```
to do a rocking scan or other sort of scan. Due to the communication overhead, sometimes it may be useful to instead run the rocking curve binning on just the psana node in the way above (without the last --num_events parameter) if the scan is very short.

The above code is best to use if you need critically fast processing time during an experiment. In post beamtime processing, I suggest using only 6 nodes without the exclusive or the ntasks-per-node option above. These last two options force the process to be able to use the entire memory contents of a node.

If you have questions, you can reach out to krapivin at [inset commonly used name of Leland Stanford Jr. University].edu

## Future plans
1. More command line options
2. Additional batch scripts to quickly run the cube code with certain parameters set
3. Removal of all hard code event codes.
4. Rebasing the entire codebase so that a particular sort of cube is actually specified in bash parameters.
5. Documentation of other code in this directory as well as command line parameters.
6. Pattern based cubing; ie. support only cubing ever 3rd shot or some other patterns.
7. Automatic generation of debug parameters like i0 v. ROI (automatically found) and other things.
