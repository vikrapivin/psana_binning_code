# script to submit for processing automatically from run table
import subprocess




sbatch -n 6 --exclusive -p psanaq process_cube.sh 351 -10 25