import psana
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import pylab as pyl
#pyl.interactive(True)
from socket import gethostname
import os
from datetime import datetime
import time
import h5py

import argparse
from bin_events_mod import binEvents
# from events_stats import collect_stats
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
savepath = '/reg/d/psdm/xpp/' + expname + '/results/krapivin/runs/r%s_laserOff.h5'%run_num
ipmlower = 1000.
ipmupper = 60000.
ttamplower = 0.01
laseroffevr = 91
laseronevr = 90
thresholdVal = 265

# get some scan parameters
if rank == 0:
  ds_get_evt = psana.DataSource('exp=' + expname + ':run=%s:smd'%run_num)
  epics = ds_get_evt.env().epicsStore()
  configStore = ds_get_evt.env().configStore()
  for nevent, ev in enumerate(ds_get_evt.events()):
    if nevent > 10:
      num_bins = epics.getPV('XPP:SCAN:NSTEPS').data()[0]
      range_lower = epics.getPV('XPP:SCAN:MIN00').data()[0]
      range_upper = epics.getPV('XPP:SCAN:MAX00').data()[0]
      scan_motor_name = configStore.get(psana.ControlData.Config).pvControls()[0].name()
      break
else:
  num_bins = None
  range_lower = None
  range_upper = None
  scan_motor_name = None
num_bins = comm.bcast(num_bins, root=0)
range_lower = comm.bcast(range_lower, root=0)
range_upper = comm.bcast(range_upper, root=0)
scan_motor_name = comm.bcast(scan_motor_name, root=0)


#num_bins = 21
#weightby = False
#control_name = None
#################################################

ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd'%run_num)
# psana.DataSource('dir=/reg/d/ffb/xpp/xpph6615/xtc/:exp=xpph6615:run=%s:idx'%run)

epics = ds.env().epicsStore()
configStore = ds.env().configStore()

evr2skip = laseronevr
total_events = 0

# ipm_src = psana.Source('BldInfo(XppSb3_Ipm)')
ipm2_src = psana.Source('BldInfo(XPP-SB2-BMMON)')
ipm3_src = psana.Source('BldInfo(XPP-SB3-BMMON)')
evr_src = psana.Source("DetInfo(NoDetector.0:Evr.0)")
evrDet   = psana.Detector('evr0',ds.env())
cspadDet = psana.Detector('jungfrau1M',ds.env()) 
encoder_src = psana.Source('DetInfo(XppEndstation.0:USDUSB.0)')

bin_ipm2 = binEvents.init_from_array(np.linspace(range_lower,range_upper,num_bins))
bin_ipm3 = binEvents.init_from_array(np.linspace(range_lower,range_upper,num_bins))
bin_cspad_sum = binEvents.init_from_array(np.linspace(range_lower,range_upper, num_bins))

def mpi_message(msg):
  if rank == 0:
    print(msg)

#def smdgen(evts):
#  for nevent, ev in enumerate(evts):
#    if nevent%size == rank: yield ev

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

    evr=ev.get(psana.EvrData.DataV4, evr_src)
    if evr is None:
      print("*** no evr, skipping.")
      continue
    evr_list = [x.eventCode() for x in evr.fifoEvents()]
    if evr_list.count(evr2skip)>0:
      skipctr += 1
      continue
    # yield ev

    # TTvalue = epics.getPV('XPP:TIMETOOL:FLTPOS_PS').data()[0]
    # TTampl = epics.getPV('XPP:TIMETOOL:AMPL').data()[0]
    # if TTampl < ttamplower:
    #   skipctr += 1
    #   # print "*** skipping TTvalue."
    #   continue
    yield ev


def get_control_value():
  control = configStore.get(psana.ControlData.Config).pvControls()
  # precondition: len(control) == 1
  if control[0].name() == "lxt_ttc" or control[0].name() == "lxt" or control[0].name() == "lxt_new" or control[0].name() == "xpp_lxt_fast":
    raise("This is a delay scan. Handle this error or comment it if you are sure you would like to continue.")
    return 1e12*control[0].value()
  else:
    return control[0].value()


## MPI reduce (only for arrays):
def mpi_reduce_arrays(*args):
  out = []
  for x in args:
    if x is not None:
      tmp = np.zeros(x.shape, dtype = x.dtype)
      comm.Reduce(x, tmp)
    else: 
      tmp = x
    out.append(tmp)
  return out

print(num_bins)


def get_config_string():
  msg = """
  expname = {}
  savepath = {}
  ipmlower = {}
  ipmupper = {}
  laseroffevr = {}
  laseronevr = {}
  range_lower = {}
  range_upper = {}
  scan_motor_name = {}
  num_bins = {}
  thresholdVal = {}
    """
  s = msg.format( expname,
                  savepath,
                  ipmlower,
                  ipmupper,
                  laseroffevr,
                  laseronevr,
                  range_lower,
                  range_upper,
                  scan_motor_name,
                  num_bins,
                  thresholdVal)
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

for nevent, ev in enumerate(filter_events(ds.events() )):
  total_events += 1
  # maybe this goes after ipm's to capture laseroffs at the end:
  scan_val = get_control_value()
  if scan_val is None:
    continue

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
    cspad_data[cspad_data <= thresholdVal] = 0
    bin_cspad_sum.update_bins(scan_val, cspad_data)
    bin_ipm2.update_bins(scan_val, ipm2intens)
    bin_ipm3.update_bins(scan_val, ipm3intens)

  if (nevent%500==0):
    # print("# events: ", nevent)
    mpi_message("number of events: %s"% nevent)
    mpi_message("stage position = %s"% scan_val)

bin_cspad_mean, bin_cspad_sum_count, bin_ipm2_mean, bin_ipm2_cts, bin_ipm3_mean, bin_ipm3_cts = mpi_reduce_arrays(bin_cspad_sum._img, 
                    bin_cspad_sum._bin_count, 
                    bin_ipm2._img, 
                    bin_ipm2._bin_count,
                    bin_ipm3._img, 
                    bin_ipm3._bin_count)

scan_vals = bin_cspad_sum.bin_edges()

if rank==0:
  saveh5(savepath,  cspad = bin_cspad_mean, 
                    nEntries = bin_cspad_sum_count, 
                    ipm2_sum = bin_ipm2_mean, 
                    ipm3_sum = bin_ipm3_mean, 
                    bins = scan_vals,
                    config = get_config_string())
  print("****** Done. ")
