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
parser.add_argument("--num_events_start", help="event index to begin processing at.", type=int, default=0)
parser.add_argument("--laser_off", help="whether or not to add laser offs to the cube", type=int, default=0)
parser.add_argument("--pull_from_ffb", help='pull xtc files from ffb, (1 if true, 0 if false, default false)', type=int, default=0)
parser.add_argument("--time_low", help="specify time start", type=float, default=0.0)
parser.add_argument("--time_high", help="specify time end", type=float, default=0.0)
parser.add_argument("--exp_name", help="specify the name of the experiment", type=str, default='xpplw8419')
parser.add_argument("--ipm_lower", help="lower value of the IPM to filter by", type=float, default=0.0)
parser.add_argument("--ipm_higher", help="higher value of the IPM to filter by", type=float, default=1.0e10)
parser.add_argument("--ipm3_lower", help="lower value of the IPM 3 to filter by", type=float, default=0.0)
parser.add_argument("--ipm3_higher", help="higher value of the IPM 3 to filter by", type=float, default=1.0e10)
parser.add_argument("--tt_amp", help="time tool amplitude filter (ignore shots below this amplitude).", type=float, default=0.001)
parser.add_argument("--bin_number", help="total number of bins of make", type=int, default=400)
parser.add_argument("--detector_threshold", help="lower value to threshold the detector (ie. values below this are set to 0)", type=float, default=0.0)
parser.add_argument("--detector_threshold_high", help="higher value to threshold the detector (ie. values above this are set to 0)", type=float, default=1.0e10)
parser.add_argument("--tt_calibration", help="The slope of the time tool calibration output, this is in units of picoseconds per pixel", type=float, default=-0.00156375)
parser.add_argument("--tt_offset", help="The offset of the time tool calibration line. Not necessary, but then the times will be shifted by the average of the time tool.", type=float, default=0.88457732)
parser.add_argument("--use_fitted_tt", help="Set to a value greater than 0 to use the fitted timetool pixel position, provided by the facility. By default the facility fit will not be used.", type=int, default=0)
parser.add_argument("--custom_calibration_directory", help="Provide a custom directory for the pedestal and other psana settings", type=str, default='')
parser.add_argument("--save_directory", help="Provide a directory (from the root of the experimental folder) to save the cube in", type=str, default='/results/krapivin/runs/')
parser.add_argument("--append_file_name", help="Append a string to the filename in the save directory", type=str, default='')


args = parser.parse_args()
run_num = args.run
num_events_limit = args.num_events
num_events_start = args.num_events_start
if size > 1: # ugly hack to do event skipping in MPI mode. As long as you are doing approximate numbers, this is approximately correct.
  num_events_limit = int(num_events_limit/size)
  num_events_start = int(num_events_start/size)
process_laser_off = args.laser_off
pull_from_ffb = args.pull_from_ffb
expname = args.exp_name
ipmlower = args.ipm_lower
ipmupper = args.ipm_higher
ipm3lower = args.ipm3_lower
ipm3upper = args.ipm3_higher
ttamplower = args.tt_amp
num_bins = args.bin_number
thresholdVal = args.detector_threshold
thresholdVal_max = args.detector_threshold_high
TIME_TOOL_CALIB = args.tt_calibration
TIME_TOOL_OFFSET = args.tt_offset
saveDirectory = args.save_directory
appendFilename = "".join(x for x in (args.append_file_name) if x.isalnum()) # make sure string appending is a valid thing to append

use_fitted_tt = True if args.use_fitted_tt > 0 else False

if pull_from_ffb == 0:
  # savepath = '/reg/d/psdm/xpp/' + expname + '/results/krapivin/runs/r%s.h5'%run_num
  savepath = '/reg/d/psdm/xpp/' + expname + saveDirectory +('r%s'%run_num) + appendFilename + '.h5'
else:
  savepath = '/cds/data/drpsrcf/xpp/'+ expname + '/scratch/krapivin/runs/r%s.h5'%run_num


laseroffevr = 91
laseronevr = 90
delay_range_lower = args.time_low
delay_range_upper = args.time_high # in units of ps

# count how many events are skipped
skipctr = 0
count = 0


#ds = psana.DataSource('exp=' + expname + ':run=%s:smd'%run_num)
# ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd'%run_num)
# psana.DataSource('dir=/reg/d/ffb/xpp/xpph6615/xtc/:exp=xpph6615:run=%s:idx'%run)
if args.custom_calibration_directory != '':
  print('Using a custom calibration directory: ' + args.custom_calibration_directory)
  psana.setOption('psana.calib-dir',args.custom_calibration_directory)

if pull_from_ffb == 0:
  ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd'%run_num)
else:
  ds = psana.MPIDataSource('exp=' + expname + ':run=%s:smd:dir=/cds/data/drpsrcf/xpp/'%run_num + expname + '/xtc')

epics = ds.env().epicsStore()
configStore = ds.env().configStore()

evr2skip = laseroffevr
total_events = 0

ipm2_src = psana.Source('BldInfo(XPP-SB2-BMMON)')
ipm3_src = psana.Source('BldInfo(XPP-SB3-BMMON)')
evr_src = psana.Source("NoDetector.0:Evr.0")
cspadDet =  psana.Detector('jungfrau1M',ds.env()) 
encoder_src = psana.Detector('XPP-USB-ENCODER-02')

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
  global skipctr
  global count
  # count2 = 0
  for ev in evts:
    # print(skipctr)
    count += 1
    if count < num_events_start:
      skipctr += 1
      continue
    if count > num_events_limit:
      print(f'processed: {count}, events.')
      print(f'skipped: {skipctr}, events.')
      break

    ## Filter on IPM2
    if ipmlower > 0:
      evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1, ipm2_src )
      if evt_intensity is not None:
        intens_ = evt_intensity.TotalIntensity()
        if intens_ < ipmlower or intens_ > ipmupper:
          skipctr += 1
          continue

    ## Filter on IPM3
    if ipm3lower > 0:
      evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1, ipm3_src )
      if evt_intensity is not None:
        intens_ = evt_intensity.TotalIntensity()
        if intens_ < ipm3lower or intens_ > ipm3upper:
          skipctr += 1
          continue

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
  timetool_delay = epics.getPV('XPP:TIMETOOL:FLTPOS').data()[0] 
  timetool_ampl  = epics.getPV('XPP:TIMETOOL:AMPL').data()[0]
  return timetool_delay, timetool_ampl

def get_timetool_facility_converted_values():
  timetool_delay = epics.getPV('XPP:TIMETOOL:FLTPOS_PS').data()[0] 
  timetool_ampl  = epics.getPV('XPP:TIMETOOL:AMPL').data()[0]
  return timetool_delay, timetool_ampl

def get_corrected_delay(ev):
  nominal_delay = get_encoder_value(ev)
  if nominal_delay is None:
    return None
  if use_fitted_tt:
    timetool_delay, timetool_ampl = get_timetool_values()
    return nominal_delay + timetool_delay
  else:
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
  run_num = {}
  savepath = {}
  ipmlower = {}
  ipmupper = {}
  ipm3lower = {}
  ipm3upper = {}
  ttamplower = {}
  TIME_TOOL_CALIB = {}
  TIME_TOOL_OFFSET = {}
  laseroffevr = {}
  laseronevr = {}
  delay_range_lower = {}
  delay_range_upper = {}
  num_bins = {}
  thresholdVal = {}
  thresholdVal_max = {}
  process_laser_off = {}
  pull_from_ffb = {}
  use_fitted_tt = {}
  num_events_start = {}
  num_events_limit = {}
  """
  s = msg.format( expname,
                  run_num,
                  savepath,
                  ipmlower,
                  ipmupper,
                  ipm3lower,
                  ipm3upper,
                  ttamplower,
                  TIME_TOOL_CALIB,
                  TIME_TOOL_OFFSET,
                  laseroffevr,
                  laseronevr,
                  delay_range_lower,
                  delay_range_upper,
                  num_bins,
                  thresholdVal,
                  thresholdVal_max,
                  process_laser_off,
                  pull_from_ffb,
                  use_fitted_tt,
                  num_events_start,
                  num_events_limit)
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

ipm3dataexists = False
ipm2dataexists = False

# for nevent, ev in enumerate(filter_events(smdgen(ds.events() ))):
for nevent, ev in enumerate(filter_events(ds.events() )):
  total_events += 1

  evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1,ipm2_src )
  if evt_intensity is None:
    if ipmlower < 0 and ipm2dataexists == False:
      ipm2intens = 0
    else:
      print("*** missing Ipm2 data. Skipping event...")
      continue
  else:
    ipm2intens = evt_intensity.TotalIntensity()
    ipm2dataexists = True
  
  evt_intensity=ev.get(psana.Bld.BldDataBeamMonitorV1, ipm3_src )
  if evt_intensity is None:
    if ipm3lower < 0 and ipm3dataexists == False:
      ipm3intens = 0
    else:
      print("*** missing Ipm3 data. Skipping event...")
      continue
  else:
    ipm3intens = evt_intensity.TotalIntensity()
    ipm3dataexists = True
  
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

#delays = bin_cspad_sum.bin_edges()
#delays = bin_cspad_sum.bin_centers()
delays = np.linspace(delay_range_lower, delay_range_upper, num_bins,endpoint=True)
mpi_message(bin_cspad_sum_count)
# print(f'\ncur thread: Filtered {skipctr} events out of a total of {count} events')
count_f = 0
skipctr_f = 0
count_f = comm.reduce(count)
skipctr_f = comm.reduce(skipctr)
mpi_message(f'\nFiltered {skipctr_f} events out of a total of {count_f} events')

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


