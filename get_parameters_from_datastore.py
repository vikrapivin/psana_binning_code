import psana
import numpy as np
import os
import h5py
from bin_events_mod import binEvents
expname='xpplw8419'
run_num=223
ds_get_evt = psana.DataSource('exp=' + expname + ':run=%s:smd'%run_num)
epics = ds_get_evt.env().epicsStore()
configStore = ds_get_evt.env().configStore()



ipm2_src = psana.Source('BldInfo(XPP-SB2-BMMON)')
ipm3_src = psana.Source('BldInfo(XPP-SB3-BMMON)')
evr_src = psana.Source("DetInfo(NoDetector.0:Evr.0)")
cspadDet = psana.Detector('jungfrau1M',ds_get_evt.env()) 
encoder_src = psana.Detector('XPP-USB-ENCODER-02')


eventIterator = ds_get_evt.events()
curEvent = next(eventIterator)
    # for name in epics.names():
    #     try:
    #             print(f''%epics.getPV(name).data()[0])
    #     except:
    #             pass


print('The scanned variable is: {0}.'.format(epics.value('XPP:SCAN:SCANVAR00')))


getNames = ['XPP:SCAN:NSTEPS','lxt','lxt_fast', 'theta','swivel_x','swivel_z','sam_x','sam_y','sam_z','THz_v','robot_x','robot_y','robot_z','cryo_temp_sample','cryo_temp_base']
for name in getNames:
    try:
        print(name + ': ' + f'{epics.getPV(name).data()[0]}')
    except:
        pass

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

# def get_corrected_delay(ev):
#   nominal_delay = get_encoder_value(ev)
#   if nominal_delay is None:
#     return None
#   timetool_delay, timetool_ampl = get_timetool_values()
#   return nominal_delay + TIME_TOOL_CALIB*timetool_delay + TIME_TOOL_OFFSET


def get_nominal_delay(ev):
  nominal_delay = get_encoder_value(ev)
  if nominal_delay is None:
    return None
  return nominal_delay

# previousDelay = float("NAN")
# previousDelay2 = previousDelay
# breakCondition1 = False
# breakCondition2 = False
# keepSearchingOnce = True
# eventNumToKeepSearching = -1
# for nevent, ev in enumerate(ds_get_evt.events()):
#     delay = get_nominal_delay(ev)
#     if nevent == 0:
#         minDelay = delay
#         maxDelay = delay
#     else:
#         if minDelay > delay:
#             minDelay = delay
#         elif maxDelay < delay:
#             maxDelay = delay
#         else:
#             if previousDelay > delay and previousDelay2 < previousDelay:
#                 breakCondition1 = True
#             if previousDelay < delay and previousDelay2 > previousDelay:
#                 breakCondition2 = True
#     if breakCondition1 == True and breakCondition2 == True:
#         if keepSearchingOnce == True:
#             eventNumToKeepSearching = nevent
#             keepSearchingOnce = False
#         elif eventNumToKeepSearching + 1000 > nevent:
#             breakCondition1 = False
#             breakCondition2 = False
#         else:
#             break
#     # print(f'{delay} is the delay, is the previous delay {previousDelay} and the second previous delay is {previousDelay2}. {breakCondition1} is the first condition and {breakCondition2} is the second condition.')
#     previousDelay2 = previousDelay
#     previousDelay = delay
#     # evr=ev.get(psana.EvrData.DataV4, evr_src)
#     # if evr is None:
#     #     pass
#     # else:
#     #     evr_list = [x.eventCode() for x in evr.fifoEvents()]
#     # print(evr_list)
#     # print(delay

# print(f'The minimum delay is {minDelay}. The maximum delay is {maxDelay}. Number of events considered is {nevent}')


# values = []
# names = []
# for nevent, ev in enumerate(ds_get_evt.events()):
# for name in epics.names():
#   try:
#     # print(f'{epics.getPV(name).data()[0]}')
#     values.append(epics.getPV(name).data()[0])
#     names.append(name)
#   except:
#     pass

# epics.getPV('XPP:TIMETOOL:AMPLNXT').data()



# for nevent, ev in enumerate(ds_get_evt.events()):
#   evr=curEvent.get(psana.EvrData.DataV4, evr_src)
#   if evr is None:
#       pass
#   else:
#       evr_list = [x.eventCode() for x in evr.fifoEvents()]
#       for evr_code in evr_list:
#         if evr_code == 162:
#           print(evr_list)

# print(evr_list)
