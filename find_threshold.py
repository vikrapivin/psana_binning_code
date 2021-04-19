import psana
import numpy as np
import matplotlib.pyplot as plt
# matplotlib.use('Agg')
import os
# from datetime import datetime
# import time
import h5py
import argparse
from bin_events_mod import binEvents


parser = argparse.ArgumentParser()
parser.add_argument("run", help="run number or range, e.g. 100,109-113.", type=str)
parser.add_argument("--event_number", help="which event number to start for threshold", type=int, default=0)
parser.add_argument("--ymin", help="which event number to start for threshold", type=int, default=0)
parser.add_argument("--ymax", help="which event number to start for threshold", type=int, default=100)
parser.add_argument("--xmin", help="which event number to start for threshold", type=int, default=0)
parser.add_argument("--xmax", help="which event number to start for threshold", type=int, default=1000)
parser.add_argument("--bins", help="which event number to start for threshold", type=int, default=40)
parser.add_argument("--savefig", help="which event number to start for threshold", type=int, default=1)
args = parser.parse_args()
run_num = args.run
start_event_num = args.event_number
ymin = args.ymin
ymax = args.ymax
xmin = args.xmin
xmax = args.xmax
bins_count = args.bins
saveFile = args.savefig
########## set parameters here: #################
expname = 'xpplv2818'
savepath = '/reg/d/psdm/xpp/' + expname + '/results/krapivin/runs/r%s.h5'%run_num
ipmlower = 1000.0
ipmupper = 60000.0
laseroffevr = 91
laseronevr = 90
########## end parameters #################

ds = psana.DataSource('exp=' + expname + ':run=%s:smd'%run_num)
epics = ds.env().epicsStore()
configStore = ds.env().configStore()

# ipm_src = psana.Source('BldInfo(XppSb3_Ipm)')
ipm2_src = psana.Source('BldInfo(XPP-SB2-BMMON)')
ipm3_src = psana.Source('BldInfo(XPP-SB3-BMMON)')
evr_src = psana.Source("DetInfo(NoDetector.0:Evr.0)")
cspadDet = psana.Detector('jungfrau1M',ds.env()) 


""" Filters by skipping events that are outside of the set of events with the desired properties. Comment out what you do not want.
"""
def filter_events(evts):
  skipctr = 0
  count = 0
  for ev in evts:
    count += 1

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
    if evr_list.count(laseronevr)>0:
      skipctr += 1
      continue
    yield ev


def get_config_string(event_num_in):
  msg = """
  expname = {}
  savepath = {}
  ipmlower = {}
  ipmupper = {}
  eventNumber {}
    """
  s = msg.format( expname,
                  savepath,
                  ipmlower,
                  ipmupper,
                  event_num_in)
  return s

total_events = 0

for nevent, ev in enumerate(filter_events(ds.events() )):
  total_events += 1
  if total_events < start_event_num:
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
    print(get_config_string(nevent))
    print(ipm2intens)
    print(ipm3intens)
    plt.figure()
    plt.hist(cspad_data.flatten(),range=(xmin,xmax),bins=bins_count)
    plt.ylim(ymin,ymax)
    print(saveFile)
    if saveFile < 1:
      plt.show()
    else:
      plt.savefig('hist.png')
      break
    
