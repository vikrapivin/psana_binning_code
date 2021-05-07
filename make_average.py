import psana
import numpy as np
from socket import gethostname
import os
from datetime import datetime
import time
import h5py
import argparse
from bin_events_mod import binEvents


# from mpi4py import MPI
# comm = MPI.COMM_WORLD 
# rank = comm.Get_rank()
# size = comm.Get_size()

cspad_avg = None
ipm2_sum_avg = None
ipm3_sum_avg = None
nEntries_avg = None
bins = None
config = None


########## set parameters here: #################
expname = 'xpplv2818'
filename_prepend = '/reg/d/psdm/xpp/' + expname + '/results/krapivin/runs/'
################################################


runs_to_average = np.array([72, 73, 74, 75, 76, 77, 78, 79])
fileNameAppend = ''
for run_num in runs_to_average:
  loadpath = filename_prepend + 'r%s.h5'%run_num
  fileNameAppend += str(run_num)
  #open file
  curRun = h5py.File(loadpath, "r")
  if cspad_avg is None:
    cspad_avg = curRun['cspad'][()]
    ipm2_sum_avg = curRun['ipm2_sum'][()]
    ipm3_sum_avg = curRun['ipm3_sum'][()]
    nEntries_avg = curRun['nEntries'][()]
    bins = curRun['bins'][()]
    config = "".join(map(chr, curRun['config'][()])) + 'Average of Runs: ' + str(runs_to_average)
    print(config)
  else:
    cspad_avg += curRun['cspad'][()]
    ipm2_sum_avg += curRun['ipm2_sum'][()]
    ipm3_sum_avg += curRun['ipm3_sum'][()]
    nEntries_avg += curRun['nEntries'][()]
  


## saving variables to HDF5
def saveh5(fullpath, **kwargs):
  """Saves a bunch of variables with names passed through kwargs dict to an HDF5 file.
  Utility h5dump comes handy when you need to look quickly into the HDF5 contents from the command line.
  """
  fulldir = os.path.dirname(fullpath)
  if (not os.path.exists(fulldir)):
    if not fulldir == '':
      try:
        print("creating directory " + fulldir)
        os.makedirs(fulldir)
      except:
        print("Could not create directory ", fulldir, "\n giving up...")
        raise
  print("Saving configuration and the following variables: ")
  for k in kwargs.keys():
    print(k)
  print("To file : %s" % fullpath)
  fid = h5py.File(fullpath, "w")
  for k, v in kwargs.items():
    if v is None:
      print("***** Warning: found empty output for " + str(k ))
    else:
      fid[k] = v
  fid.close()

savepath = filename_prepend + 'r%s_avg.h5'%fileNameAppend
saveh5(savepath,  cspad = cspad_avg, 
                  nEntries = nEntries_avg, 
                  ipm2_sum = ipm2_sum_avg, 
                  ipm3_sum = ipm3_sum_avg, 
                  bins = bins,
                  config = config)
