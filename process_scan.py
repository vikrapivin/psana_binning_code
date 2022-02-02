import psana
import numpy as np
import os
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
parser.add_argument("--exp_name", help="specify the name of the experiment", type=str, default='xpplw8419')
parser.add_argument("--background", help="name of background file, leave empty if none", type=str, default='')
parser.add_argument("--relative_scan", help="treat scan as relative scan (1 if true, 0 if false, -1 if the scan is reversed direction, default false)", type=int, default=0)
parser.add_argument("--pull_from_ffb", help='pull xtc files from ffb, (1 if true, 0 if false, default false)', type=int, default=0)
parser.add_argument("--laser_off", help="whether or not to add a laser off cube to the cube", type=int, default=0)
parser.add_argument("--ignore_no_optical_laser", help="Should this code explicitly check that each laser on event has a laser on code. This is typically not what you want as you want data from when we have optical laser off.", type=int, default=0)
parser.add_argument("--ipm_lower", help="lower value of the IPM to filter by (-1 to disable filter)", type=float, default=0.0)
parser.add_argument("--ipm_higher", help="higher value of the IPM to filter by", type=float, default=1.0e10)
parser.add_argument("--ipm3_lower", help="lower value of the IPM 3 to filter by (-1 to disable filter)", type=float, default=0.0)
parser.add_argument("--ipm3_higher", help="higher value of the IPM 3 to filter by", type=float, default=1.0e10)
parser.add_argument("--detector_threshold", help="lower value to threshold the detector (ie. values below this are set to 0)", type=float, default=0.0)
parser.add_argument("--detector_threshold_high", help="higher value to threshold the detector (ie. values above this are set to 0)", type=float, default=1.0e10)
parser.add_argument("--custom_calibration_directory", help="Provide a custom directory for the pedestal and other psana settings", type=str, default='')
parser.add_argument("--save_directory", help="Provide a directory (from the root of the experimental folder) to save the cube in", type=str, default='/results/krapivin/runs/')
parser.add_argument("--append_file_name", help="Append a string to the filename in the save directory", type=str, default='')


# TODO: add parameters that corrrespond to detectors, psana source, etc.

args = parser.parse_args()
run_num = args.run
num_events_limit = args.num_events
background_file = args.background
is_relative_scan = args.relative_scan
pull_from_ffb = args.pull_from_ffb
process_laser_off = args.laser_off
ignore_no_optical_laser = args.ignore_no_optical_laser
expname = args.exp_name
ipmlower = args.ipm_lower
ipmupper = args.ipm_higher
ipm3lower = args.ipm3_lower
ipm3upper = args.ipm3_higher
thresholdVal = args.detector_threshold
thresholdVal_max = args.detector_threshold_high
saveDirectory = args.save_directory
appendFilename = "".join(x for x in (args.append_file_name) if x.isalnum()) # make sure string appending is a valid thing to append

########## set parameters here: #################
if pull_from_ffb == 0:
  savepath = '/reg/d/psdm/xpp/' + expname + saveDirectory +('r%s'%run_num) + appendFilename + '.h5'
else:
  # fixME
  savepath = '/cds/data/drpsrcf/xpp/'+ expname + '/scratch/krapivin/runs/r%s.h5'%run_num

laseroffevr = 91
laseronevr = 90

skipctr = 0
count = 0

# get some scan parameters
rel_offset = 0
scanMotorConfigStoreOffSet = 0
if rank == 0:
  if pull_from_ffb == 0:
    ds_get_evt = psana.DataSource('exp=' + expname + ':run=%s:smd'%run_num)
  else:
    ds_get_evt = psana.DataSource('exp=' + expname + ':run=%s:smd:dir=/cds/data/drpsrcf/xpp/'%run_num + expname + '/xtc')
  epics = ds_get_evt.env().epicsStore()
  configStore = ds_get_evt.env().configStore()
  for nevent, ev in enumerate(ds_get_evt.events()):
    if nevent > 10:
      num_bins = epics.getPV('XPP:SCAN:NSTEPS').data()[0]
      if num_bins is None:
        print('num_bins did not return. Assuming 5 bins.')
        num_bins = 5
      range_lower = epics.getPV('XPP:SCAN:MIN00').data()[0]
      range_upper = epics.getPV('XPP:SCAN:MAX00').data()[0]
      # the thing below seems to be more reliable than getting the config store as multiple scan motors can be stored in the config store
      scan_motor_name = epics.value('XPP:SCAN:SCANVAR00') #configStore.get(psana.ControlData.Config).pvControls()[0].name() 
      # print(scan_motor_name + ' is scan motor.')
      control = configStore.get(psana.ControlData.Config).pvControls()
      lengthOfSavedControls = len(control)
      # find the correct control parameter to use, assuming we only have 1 scanvar, TODO to handle if more than 1 scan var
      foundVal = False
      for ii in range(0,lengthOfSavedControls):
        # print(control[ii].name() + ' is ii motor.')
        if(control[ii].name() == scan_motor_name):
          scanMotorConfigStoreOffSet = ii
          foundVal = True
          break
      if foundVal == False:
        print('Warning: was not able to find the scan motor config store offset. Assuming it is 0.')
      if is_relative_scan >0:
        rel_offset = configStore.get(psana.ControlData.Config).pvControls()[scanMotorConfigStoreOffSet].value() - range_lower
      if is_relative_scan <0:
        rel_offset = configStore.get(psana.ControlData.Config).pvControls()[scanMotorConfigStoreOffSet].value() - range_upper
      break
  # get background if specified
  if background_file != '':
    bkgrFile = h5py.File(background_file,'r')
    avgBkgr = bkgrFile['avgBkgr'][()].astype('float64')
    cspad_data_shape = avgBkgr.shape
    bkgrFile.close()
else:
  cspad_data_shape = None
  num_bins = None
  range_lower = None
  range_upper = None
  scan_motor_name = None
num_bins = comm.bcast(num_bins, root=0)
range_lower = comm.bcast(range_lower, root=0)
range_upper = comm.bcast(range_upper, root=0)
scan_motor_name = comm.bcast(scan_motor_name, root=0)
scanMotorConfigStoreOffSet = comm.bcast(scanMotorConfigStoreOffSet, root=0)
if background_file != '':
  comm.Barrier()
  cspad_data_shape = comm.bcast(cspad_data_shape, root=0)
  if rank != 0:
    avgBkgr = np.zeros(cspad_data_shape, dtype='float64')
  comm.Barrier()
  comm.Bcast(avgBkgr, root=0)
if is_relative_scan !=0:
  rel_offset = comm.bcast(rel_offset, root=0)
  range_lower = rel_offset + range_lower
  range_upper = rel_offset + range_upper


#################################################


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
ipm3_src = psana.Source('BldInfo(XPP-SB3-BMMON)')
evr_src = psana.Source("DetInfo(NoDetector.0:Evr.0)")
evrDet   = psana.Detector('evr0',ds.env())
cspadDet = psana.Detector('jungfrau1M',ds.env()) 
encoder_src = psana.Detector('XPP-USB-ENCODER-02')

bin_ipm2 = binEvents.init_from_array(np.linspace(range_lower,range_upper,num_bins))
bin_ipm3 = binEvents.init_from_array(np.linspace(range_lower,range_upper,num_bins))
bin_cspad_sum = binEvents.init_from_array(np.linspace(range_lower,range_upper, num_bins))
if process_laser_off !=0:
  bin_cspad_sum_off = binEvents.init_from_array(np.linspace(range_lower,range_upper,num_bins))#(range_lower, range_upper, num_bins)
  bin_ipm2_off = binEvents.init_from_array(np.linspace(range_lower,range_upper,num_bins))#(range_lower, range_upper, num_bins)
  bin_ipm3_off = binEvents.init_from_array(np.linspace(range_lower,range_upper,num_bins))#(range_lower, range_upper, num_bins)


def mpi_message(msg):
  if rank == 0:
    print(msg)

#def smdgen(evts):
#  for nevent, ev in enumerate(evts):
#    if nevent%size == rank: yield ev

""" Filters by skipping events that are outside of the desired range. Comment out what you do not want.
"""
def filter_events(evts):
  global skipctr
  global count
  for ev in evts:
    count += 1
    if count > num_events_limit:
      print("processed: ", count, "events.")
      print("skipped: ", skipctr, "events.")
      break


    ## Filter on IPM2
    if ipmlower > 0:
      evt_intensity=ev.get(psana.Lusi.IpmFexV1, ipm2_src )
      if evt_intensity is not None:
        intens_ = evt_intensity.sum()
        if intens_ < ipmlower or intens_ > ipmupper:
          skipctr += 1
          continue

    ## Filter on IPM3
    if ipm3lower > 0:
      evt_intensity=ev.get(psana.Lusi.IpmFexV1, ipm3_src )
      if evt_intensity is not None:
        intens_ = evt_intensity.sum()
        if intens_ < ipm3lower or intens_ > ipm3upper:
          skipctr += 1
          continue

    evr=ev.get(psana.EvrData.DataV4, evr_src)
    if evr is None:
      print("*** no evr, skipping.")
      continue
    if process_laser_off !=0:
      pass
    else:
      # evr_list = [x.eventCode() for x in evr.fifoEvents()]
      # print(evr_list)
      # if evr_list.count(evr2skip)>0:
      #   skipctr += 1
      #   continue
      evr_list = [x.eventCode() for x in evr.fifoEvents()]
      if evr_list.count(laseroffevr)>0:
        skipctr += 1
        continue
      elif ignore_no_optical_laser > 0:
        if evr_list.count(laseronevr)>0:
          pass
        else:
          skipctr += 1
          continue
      else:
        pass
    yield ev


def get_control_value(ev):
  control = configStore.get(psana.ControlData.Config).pvControls()
  # precondition: len(control) == 1
  if control[0].name() == "lxt_ttc" or control[0].name() == "lxt" or control[0].name() == "lxt_new" or control[0].name() == "xpp_lxt_fast":
    raise("This is a delay scan. Handle this error or comment it if you are sure you would like to continue. You will likely need to specify custom ranges for this binning.")
    if encoder_src.values(ev) is None:
      return None
    enc_data = encoder_src.values(ev)[0]
    return enc_data
  else:
    return control[scanMotorConfigStoreOffSet].value()


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

#print(num_bins)


def get_config_string():
  msg = """
  expname = {}
  savepath = {}
  ipmlower = {}
  ipmupper = {}
  ipm3lower = {}
  ipm3upper = {}
  laseroffevr = {}
  laseronevr = {}
  range_lower = {}
  range_upper = {}
  scan_motor_name = {}
  num_bins = {}
  thresholdVal = {}
  thresholdVal_max = {}
  background_file = {}
  scanMotorConfigStoreOffSet = {}
  process_laser_off = {}
  ignore_no_optical_laser = {}
    """
  s = msg.format( expname,
                  savepath,
                  ipmlower,
                  ipmupper,
                  ipm3lower,
                  ipm3upper,
                  laseroffevr,
                  laseronevr,
                  range_lower,
                  range_upper,
                  scan_motor_name,
                  num_bins,
                  thresholdVal,
                  thresholdVal_max,
                  background_file,
                  scanMotorConfigStoreOffSet,
                  process_laser_off,
                  ignore_no_optical_laser)
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

ipm3dataexists = False
ipm2dataexists = False

for nevent, ev in enumerate(filter_events(ds.events() )):
  total_events += 1
  # maybe this goes after ipm's to capture laseroffs at the end:
  scan_val = get_control_value(ev)
  if scan_val is None:
    continue

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
  else:
    if background_file != '':
      cspad_data = cspad_data - avgBkgr
  if process_laser_off == 0:
    cspad_data[cspad_data <= thresholdVal] = 0
    cspad_data[cspad_data >= thresholdVal_max] = 0
    bin_cspad_sum.update_bins(scan_val, cspad_data)
    bin_ipm2.update_bins(scan_val, ipm2intens)
    bin_ipm3.update_bins(scan_val, ipm3intens)
  else:
    evr=ev.get(psana.EvrData.DataV4, evr_src)
    evr_list = [x.eventCode() for x in evr.fifoEvents()]
    # print(evr_list)
    if evr_list.count(laseroffevr)>0:
      cspad_data[cspad_data <= thresholdVal] = 0
      cspad_data[cspad_data >= thresholdVal_max] = 0
      bin_cspad_sum_off.update_bins(scan_val, cspad_data)
      bin_ipm2_off.update_bins(scan_val, ipm2intens)
      bin_ipm3_off.update_bins(scan_val, ipm3intens)
    elif ignore_no_optical_laser == 0 or evr_list.count(laseronevr)>0:
      cspad_data[cspad_data <= thresholdVal] = 0
      cspad_data[cspad_data >= thresholdVal_max] = 0
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
if process_laser_off !=0:
  bin_cspad_mean_off, bin_cspad_sum_count_off, bin_ipm2_mean_off, bin_ipm2_cts_off, bin_ipm3_mean_off, bin_ipm3_cts_off = mpi_reduce_arrays(
                      bin_cspad_sum_off._img, 
                      bin_cspad_sum_off._bin_count, 
                      bin_ipm2_off._img, 
                      bin_ipm2_off._bin_count,
                      bin_ipm3_off._img, 
                      bin_ipm3_off._bin_count)
#scan_vals = bin_cspad_sum.bin_centers()
scan_vals = np.linspace(range_lower,range_upper, num_bins)

count_f = 0
skipctr_f = 0
count_f = comm.reduce(count)
skipctr_f = comm.reduce(skipctr)
mpi_message(f'\nFiltered {skipctr_f} events out of a total of {count_f} events')

if rank==0:
  print(bin_cspad_sum_count)
  if process_laser_off == 0:
    saveh5(savepath,  cspad = bin_cspad_mean, 
                      nEntries = bin_cspad_sum_count, 
                      ipm2_sum = bin_ipm2_mean, 
                      ipm3_sum = bin_ipm3_mean, 
                      bins = scan_vals,
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
                      bins = scan_vals,
                      config = get_config_string())
  print("****** Done. ")

comm.Barrier()

