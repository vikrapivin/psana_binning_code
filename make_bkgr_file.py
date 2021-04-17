import psana
import numpy as np
import os
from datetime import datetime
import time
import h5py
import argparse
from bin_events_mod import binEvents
from mpi4py import MPI

comm = MPI.COMM_WORLD 
rank = comm.Get_rank()
size = comm.Get_size()

print(f'rank {rank} out of {size}')

parser = argparse.ArgumentParser()
parser.add_argument("run", help="run number or range, e.g. 100,109-113.", type=str)
parser.add_argument("--num_events", help="number of events to process", type=int, default=1<<31)
args = parser.parse_args()
run_num = args.run
num_events_limit = args.num_events

########## set parameters here: #################
expname = 'xpplv2818'
savepath = '/reg/d/psdm/xpp/' + expname + '/results/krapivin/runs/r%s.h5'%run_num
ipmlower = 1000.
ipmupper = 60000.
ttamplower = 0.01
#################################################


# get detector size
if rank == 0:
  ds_get_evt = psana.DataSource('exp=' + expname + ':run=%s:smd'%run_num)
  cspadDet = psana.Detector('jungfrau1M',ds_get_evt.env()) 
  cspad_data_shape = None
  for nevent, ev in enumerate(ds_get_evt.events()):
    if nevent > 10:
      cspad_data=cspadDet.image(ev)
      if cspad_data is None:
        print("*** missing cspad data. Skipping event...")
        continue
      break
  cspad_data_shape = cspad_data.shape
else:
  cspad_data_shape = None

# print(cspad_data_shape)
comm.Barrier()
cspad_data_shape = comm.bcast(cspad_data_shape, root=0)
####

ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd'%run_num)
# psana.DataSource('dir=/reg/d/ffb/xpp/xpph6615/xtc/:exp=xpph6615:run=%s:idx'%run)

epics = ds.env().epicsStore()
configStore = ds.env().configStore()
total_events = 0
cspadDet = psana.Detector('jungfrau1M',ds.env()) 


def mpi_message(msg):
  if rank == 0:
    print(msg)

""" Filters by skipping events that are outside of the desired range. Comment out what you do not want.
"""
def filter_events(evts):
  skipctr = 0
  count = 0
  for ev in evts:
    count += 1
    if count > num_events_limit:
      print("processed: ", count, "events.")
      print("skipped: ", skipctr, "events.")
      break
    yield ev


def get_config_string():
  msg = """
  expname = {}
  savepath = {}
    """
  s = msg.format( expname,
                  savepath)
  return s


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
      print("***** Warning: found empty output for " + k )
    else:
      fid[k] = v
  fid.close()


##### Main loop : 
# NOTE: if processing multiple runs, can access each run separately by iterating over ds.runs() before iterating over events.
# NOTE2: iterate over run.steps() to access steps (CalibCycles)
#  if we don't specify calib_cycles or runs we get all events from all requested runs sequentially.
# runs = ds.runs()
# for r in runs:

# event loop starts here
mpi_message("Starting event loop. ")
mpi_message("Using config: ")
mpi_message(get_config_string())

bkgrSum = np.zeros(cspad_data_shape)

for nevent, ev in enumerate(filter_events(ds.events() )):
  total_events += 1
  cspad_data=cspadDet.image(ev)
  if cspad_data is None:
    print("*** missing cspad data. Skipping event...")
    continue
  else:
    bkgrSum +=cspad_data

  if (nevent%50==0):
    # print("# events: ", nevent)
    mpi_message("number of events: %s"% nevent)
bkgrSum_final = np.zeros(bkgrSum.shape)
total_events_final = 0
comm.Reduce(bkgrSum, bkgrSum_final, root=0)
total_events_final = comm.reduce(total_events, root=0)
if rank==0:
  avgBkgr = bkgrSum_final/total_events_final
  saveh5(savepath,  cspad = bkgrSum_final, 
                    avgBkgr = avgBkgr,
                    nEntries = total_events_final, 
                    config = get_config_string())
  print("****** Done. ")
