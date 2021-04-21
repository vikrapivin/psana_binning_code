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
ipmlower = 2000.
ipmupper = 60000.
laseroffevr = 91
laseronevr = 90
thresholdVal = 9.75
thresholdVal_max = 11.0
#evr2skip = laseroffevr
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
ipm2_src = psana.Source('BldInfo(XPP-SB2-BMMON)')
ipm3_src = psana.Source('BldInfo(XPP-SB3-BMMON)')
evr_src = psana.Source("DetInfo(NoDetector.0:Evr.0)")


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
    ## Filter on IPM2
    evt_intensity=ev.get(psana.Lusi.IpmFexV1, ipm2_src )
    if evt_intensity is not None:
      intens_ = evt_intensity.sum()
      if intens_ < ipmlower or intens_ > ipmupper:
        skipctr += 1
        continue

    # ## Filter on IPM3
    # evt_intensity=ev.get(psana.Lusi.IpmFexV1, ipm3_src )
    # if evt_intensity is not None:
    #   intens_ = evt_intensity.sum()
    #   if intens_ < ipmlower or intens_ > ipmupper:
    #     skipctr += 1
    #     continue

    #evr=ev.get(psana.EvrData.DataV4, evr_src)
    #if evr is None:
    #  print("*** no evr, skipping.")
    #  continue
    #evr_list = [x.eventCode() for x in evr.fifoEvents()]
    #if evr_list.count(evr2skip)>0:
    #  skipctr += 1
    #  continue

    yield ev


def get_config_string():
  msg = """
  expname = {}
  savepath = {}
  ipmlower = {}
  ipmupper = {}
  thresholdVal = {}
  thresholdVal_max = {}
    """
  s = msg.format( expname,
                  savepath,
                  ipmlower,
                  ipmupper,
                  thresholdVal,
                  thresholdVal_max)
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

single_img = np.zeros(cspad_data_shape)
ipm2_sum = 0
ipm3_sum = 0
for nevent, ev in enumerate(filter_events(ds.events() )):
  total_events += 1

  evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1,ipm3_src )
  if evt_intensity is None:
    print("*** missing Ipm3 data. Skipping event...")
    continue
  ipm3intens = evt_intensity.TotalIntensity()

  evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1,ipm2_src )
  if evt_intensity is None:
    print("*** missing Ipm2 data. Skipping event...")
    continue
  ipm2intens = evt_intensity.TotalIntensity()

  cspad_data=cspadDet.image(ev)
  if cspad_data is None:
    print("*** missing cspad data. Skipping event...")
    continue
  else:
    cspad_data[cspad_data >= thresholdVal_max] = 0
    cspad_data[cspad_data <= thresholdVal] = 0
    single_img +=cspad_data
    ipm2_sum += ipm2intens
    ipm3_sum += ipm3intens

  if (nevent%50==0):
    # print("# events: ", nevent)
    mpi_message("number of events: %s"% nevent)
single_img_final = np.zeros(single_img.shape)
total_events_final = 0
ipm2_final = 0
ipm3_final = 0

comm.Reduce(single_img, single_img_final, root=0)
total_events_final = comm.reduce(total_events, root=0)
ipm2_final = comm.reduce(ipm2_sum, root=0)
ipm3_final = comm.reduce(ipm3_sum, root=0)
if rank==0:
  avgImg = single_img_final/total_events_final
  saveh5(savepath,  cspad = single_img_final, 
                    ipm2_sum = ipm2_final,
                    ipm3_sum = ipm3_final,
                    avgImg = avgImg,
                    nEntries = total_events_final, 
                    config = get_config_string())
  print("****** Done. ")
