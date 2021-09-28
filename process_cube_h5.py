import psana
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import pylab as pyl
pyl.interactive(True)
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
parser.add_argument("--laser_off", help="whether or not to add a laser off cube to the cube", type=int, default=0)
parser.add_argument("--pull_from_ffb", help='pull xtc files from ffb, (1 if true, 0 if false, default false)', type=int, default=0)
parser.add_argument("--time_low", help="specify time start", type=float, default=0.0)
parser.add_argument("--time_high", help="specify time end", type=float, default=0.0)
args = parser.parse_args()
run_num = args.run
num_events_limit = args.num_events
process_laser_off = args.laser_off
pull_from_ffb = args.pull_from_ffb

########## set parameters here: #################
expname = 'xpplv9818'
if pull_from_ffb == 0:
  savepath = '/reg/d/psdm/xpp/' + expname + '/results/krapivin/runs/r%s.h5'%run_num
else:
  savepath = '/cds/data/drpsrcf/xpp/'+ expname + '/scratch/krapivin/runs/r%s.h5'%run_num
#savepath = 'r%s.h5'%run_num
ipmlower = 50.0
ipmupper = 60000.0
ttamplower = 0.015
#ttfltposps_lower=
#ttfltposps_upper=

thresholdVal = 11.0
thresholdVal_max = 110000.0

TIME_TOOL_CALIB = -0.0019428
TIME_TOOL_OFFSET = 0.90017463
#end

laseroffevr = 91
laseronevr = 90
delay_range_lower = args.time_low
delay_range_upper = args.time_high # in units of ps
## use these to skip out of range delays:
delay_ignore_lower = delay_range_lower   
delay_ignore_upper = delay_range_upper 
num_bins = 400
weightby = False
encoder_conv = 0.13333333333333333e-3 # ps/click
encoder_offset= -182700
################################################

#ds = psana.DataSource('exp=' + expname + ':run=%s:smd'%run_num)
# ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd'%run_num)
# psana.DataSource('dir=/reg/d/ffb/xpp/xpph6615/xtc/:exp=xpph6615:run=%s:idx'%run)

if pull_from_ffb == 0:
  ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd'%run_num)
else:
  ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd:dir=/cds/data/drpsrcf/xpp/'%run_num + expname + '/xtc')

epics = ds.env().epicsStore()
configStore = ds.env().configStore()

evr2skip = laseroffevr
total_events = 0

# ipm_src = psana.Source('BldInfo(XppSb3_Ipm)')
ipm2_src = psana.Source('BldInfo(XPP-SB2-BMMON)')
#ipm2_src = psana.Detector('XPP-SB2-BMMON',ds.env())
ipm3_src = psana.Source('BldInfo(XPP-SB3-BMMON)')
#ipm3_src = psana.Detector('XPP-SB3-BMMON',ds.env())
evr_src = psana.Source("NoDetector.0:Evr.0")
# evrDet   = psana.Detector('evr0',ds.env())
cspadDet =  psana.Detector('jungfrau1M',ds.env()) 
#encoder_src = psana.Source('DetInfo(XppEndstation.0:USDUSB.0)')
encoder_src = psana.Detector('XPP-USB-ENCODER-02')
#encoder_src = psana.Detector()
#encoder_src = '';

bin_ipm2 = binEvents(delay_range_lower, delay_range_upper, num_bins)
bin_ipm3 = binEvents(delay_range_lower, delay_range_upper, num_bins)
bin_cspad_sum = binEvents(delay_range_lower, delay_range_upper, num_bins)
if process_laser_off !=0:
  bin_cspad_sum_off = binEvents(delay_range_lower, delay_range_upper, num_bins)
  bin_ipm2_off = binEvents(delay_range_lower, delay_range_upper, num_bins)
  bin_ipm3_off = binEvents(delay_range_lower, delay_range_upper, num_bins)


def mpi_message(msg):
  if rank == 0:
    print(msg)

# def smdgen(evts):
#   for nevent, ev in enumerate(evts):
#     if nevent%size == rank: yield ev

""" Filters by skipping events that are outside of the desired range. Comment out what you do not want.
"""
def filter_events(evts):
  skipctr = 0
  count = 0
  # count2 = 0
  for ev in evts:
    # print(skipctr)
    count += 1
    if count > num_events_limit:
      print(f'processed: {count}, events.')
      print(f'skipped: {skipctr}, events.')
      break

#    delay = get_corrected_delay(ev)
#    if delay is None or delay < delay_ignore_lower or delay > delay_ignore_upper:
#      skipctr += 1
#      continue

    ## Filter on IPM2
    evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1, ipm2_src )
    if evt_intensity is not None:
      intens_ = evt_intensity.TotalIntensity()
      # print(intens_)
      if intens_ < ipmlower or intens_ > ipmupper:
        # print('case intensity')
        skipctr += 1
        continue

    # ## Filter on IPM3
    # evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1, ipm3_src )
    # if evt_intensity is not None:
    #   intens_ = evt_intensity.TotalIntensity()
    #   if intens_ < ipmlower or intens_ > ipmupper:
    #     skipctr += 1
    #     continue

    evr=ev.get(psana.EvrData.DataV4, evr_src)
    if evr is None:
      print("*** no evr, skipping.")
      continue
    evr_list = [x.eventCode() for x in evr.fifoEvents()]
    # print(evr_list)
    if process_laser_off !=0:
      if evr_list.count(laseronevr)>0:
        # TTvalue = epics.getPV('XPP:TIMETOOL:FLTPOS').data()[0]
        TTampl = epics.getPV('XPP:TIMETOOL:AMPL').data()[0]
        if TTampl < ttamplower:
          # print(TTampl)
          # print('case TT')
          skipctr += 1
          # print "*** skipping TTvalue."
          continue
    elif evr_list.count(evr2skip)>0:
      # print(evr_list.count(evr2skip))
      skipctr += 1
      continue
    else:
      # TTvalue = epics.getPV('XPP:TIMETOOL:FLTPOS').data()[0]
      TTampl = epics.getPV('XPP:TIMETOOL:AMPL').data()[0]
      if TTampl < ttamplower:
        # print(TTampl)
        # print('case TT')
        skipctr += 1
        # print "*** skipping TTvalue."
        continue
    # count2 += 1
    # print(count2)
    yield ev

def get_encoder_value(ev):
  #enc_data = ev.get(psana.UsdUsb.DataV1, encoder_src)
  if encoder_src.values(ev) is None:
    return None
  enc_data = encoder_src.values(ev)[0]
  if enc_data is None:
    return None
  thisDat = enc_data
  return thisDat

def get_timetool_values():
  # precondition: PV.data() != None
  timetool_delay = epics.getPV('XPP:TIMETOOL:FLTPOS').data()[0] 
  timetool_ampl  = epics.getPV('XPP:TIMETOOL:AMPL').data()[0]
  return timetool_delay, timetool_ampl

def get_corrected_delay(ev):
  nominal_delay = get_encoder_value(ev)
  if nominal_delay is None:
    return None
  timetool_delay, timetool_ampl = get_timetool_values()
  return nominal_delay + TIME_TOOL_CALIB*timetool_delay + TIME_TOOL_OFFSET


def get_nominal_delay(ev):
  nominal_delay = get_encoder_value(ev)
  if nominal_delay is None:
    return None
  return nominal_delay


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


def get_config_string():
  msg = """
  expname = {}
  savepath = {}
  ipmlower = {}
  ipmupper = {}
  ttamplower = {}
  TIME_TOOL_CALIB = {}
  laseroffevr = {}
  laseronevr = {}
  delay_range_lower = {}
  delay_range_upper = {}
  delay_ignore_lower = {}
  delay_ignore_upper = {}
  num_bins = {}
  weightby = {}
  thresholdVal = {}
  thresholdVal_max = {}
  process_laser_off = {}
  """
  s = msg.format( expname,
                  savepath,
                  ipmlower,
                  ipmupper,
                  ttamplower,
                  TIME_TOOL_CALIB,
                  laseroffevr,
                  laseronevr,
                  delay_range_lower,
                  delay_range_upper,
                  delay_ignore_lower,
                  delay_ignore_upper,
                  num_bins,
                  weightby,
                  thresholdVal,
                  thresholdVal_max,
                  process_laser_off)
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
      print("***** Warning: found empty output for " + str(k ))
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

prev_delay = None

# for nevent, ev in enumerate(filter_events(smdgen(ds.events() ))):
for nevent, ev in enumerate(filter_events(ds.events() )):
  total_events += 1

  evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1,ipm2_src )
  if evt_intensity is None:
    print("*** missing Ipm2 data. Skipping event...")
    continue
  ipm2intens = evt_intensity.TotalIntensity()
  
  evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1, ipm3_src )
  if evt_intensity is None:
    print("*** missing Ipm3 data. Skipping event...")
    continue
  ipm3intens = evt_intensity.TotalIntensity()
  
  cspad_data=cspadDet.image(ev)
  if cspad_data is None:
    print("*** missing cspad data. Skipping event...")
    continue
  if process_laser_off == 0:
    delay = get_corrected_delay(ev)
    if delay is None:
      continue
    cspad_data[cspad_data <= thresholdVal] = 0
    cspad_data[cspad_data >= thresholdVal_max] = 0
    bin_cspad_sum.update_bins(delay, cspad_data)
    bin_ipm2.update_bins(delay, ipm2intens)
    bin_ipm3.update_bins(delay, ipm3intens)
  elif process_laser_off != 0:
    evr=ev.get(psana.EvrData.DataV4, evr_src)
    evr_list = [x.eventCode() for x in evr.fifoEvents()]
    # print(evr_list)
    if evr_list.count(laseroffevr)>0:
      delay = get_nominal_delay(ev)
      if delay is None:
        continue
      cspad_data[cspad_data <= thresholdVal] = 0
      cspad_data[cspad_data >= thresholdVal_max] = 0
      bin_cspad_sum_off.update_bins(delay, cspad_data)
      bin_ipm2_off.update_bins(delay, ipm2intens)
      bin_ipm3_off.update_bins(delay, ipm3intens)
    elif evr_list.count(laseronevr)>0:
      delay = get_corrected_delay(ev)
      if delay is None:
        continue
      cspad_data[cspad_data <= thresholdVal] = 0
      cspad_data[cspad_data >= thresholdVal_max] = 0
      bin_cspad_sum.update_bins(delay, cspad_data)
      bin_ipm2.update_bins(delay, ipm2intens)
      bin_ipm3.update_bins(delay, ipm3intens)
  # print(f'# events: {nevent}')
  if (total_events%100==0):
    # print(f'# events: {nevent}')
    mpi_message("# events: %s"% nevent)
# bin_count_cspad = np.copy(bin_cspad_sum._bin_count)
# comm.Reduce(bin_cspad_sum._bin_count, bin_count_cspad)
bin_cspad_mean, bin_cspad_sum_count, bin_ipm2_mean, bin_ipm2_cts, bin_ipm3_mean, bin_ipm3_cts = mpi_reduce_arrays(
                    bin_cspad_sum._img, 
                    bin_cspad_sum._bin_count, 
                    bin_ipm2._img, 
                    bin_ipm2._bin_count,
                    bin_ipm3._img, 
                    bin_ipm3._bin_count)
if process_laser_off !=0:
  bin_cspad_mean_off, bin_cspad_sum_count_off, bin_ipm2_mean_off, bin_ipm2_cts_off, bin_ipm3_mean_off, bin_ipm3_cts_off = mpi_reduce_arrays(
                      bin_cspad_sum_off._img, 
                      bin_cspad_sum_off._bin_count, 
                      bin_ipm2_off._img, 
                      bin_ipm2_off._bin_count,
                      bin_ipm3_off._img, 
                      bin_ipm3_off._bin_count)

delays = bin_cspad_sum.bin_edges()
mpi_message(bin_cspad_sum_count)
if rank==0:
  if process_laser_off == 0:
    saveh5(savepath,  cspad = bin_cspad_mean, 
                      nEntries = bin_cspad_sum_count, 
                      ipm2_sum = bin_ipm2_mean, 
                      ipm3_sum = bin_ipm3_mean, 
                      bins = delays,
                      config = get_config_string())
  elif process_laser_off != 0:
    saveh5(savepath,  cspad = bin_cspad_mean, 
                      nEntries = bin_cspad_sum_count, 
                      ipm2_sum = bin_ipm2_mean, 
                      ipm3_sum = bin_ipm3_mean, 
                      cspad_off = bin_cspad_mean_off, 
                      nEntries_off = bin_cspad_sum_count_off, 
                      ipm2_sum_off = bin_ipm2_mean_off, 
                      ipm3_sum_off = bin_ipm3_mean_off, 
                      bins = delays,
                      config = get_config_string())
  print("****** Done. ")


