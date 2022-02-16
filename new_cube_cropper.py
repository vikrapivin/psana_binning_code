import argparse
import h5py
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("run", help="run number", type=int)
parser.add_argument("--exp", help="exp name", type=str, default='xpply5120')
parser.add_argument("--expandROI", help="expand ROI by x pixels", type=int,default=0)
parser.add_argument("--roi", help="coordinates of the edges of the roi rectangle in order of: x0, x1, y0, y1", type=int, nargs=4)
args = parser.parse_args()
run = args.run
exp = args.exp
roiExp = args.expandROI
def cropCubeFileByROI(exp,run,filename,roi):
    print(f'Cropping cube to {roi} by our notation.')
    try:
        f = h5py.File(filename,'r')
    except Exception as e:
        print(f'Error opening cube file: {filename}, crop not saved for this cube.')
        return
    scan_var = f['binVar_bins'][()]
    bin_count = f['jungfrau1M_nEntries'][()]
    imgs = f['jungfrau1M_data'][:,roi[0,0]:roi[0,1],roi[1,0]:roi[1,1]]
    ROI = roi
    i0 = f['ipm2__sum'][()]
    i0_ipm3 = f['ipm3__sum'][()]
    writeFile = h5py.File(filename[0:-3]+'_cropped.h5', 'w')
    writeFile['scan_var'] = scan_var
    writeFile['bin_count'] = bin_count
    writeFile['imgs'] = imgs
    writeFile['ROI'] = ROI
    writeFile['i0'] = i0
    writeFile['i0_ipm3'] = i0_ipm3
    writeFile.close()
    return
def cropCubeFile(exp,run,filename):
    #smallDir = '/reg/data/drpsrcf/xpp/'+exp+'/scratch/hdf5/smalldata/'
    smallDir = '/cds/data/drpsrcf/xpp/'+exp+'/scratch/hdf5/smalldata/'
    smalldataName = f'xpply5120_Run{run:04}.h5'
    #print(smalldataName)
    try:
        g = h5py.File(smallDir+smalldataName,'r')
    except Exception as e:
        print(f'Error opening small data file: {smalldataName}, crop not saved for this run.')
        return
    smlROI = g['UserDataCfg/jungfrau1M/ROI_0__ROI_0_bound'][()]
    roi = np.zeros((2,2)).astype(int)
    if np.all(smlROI[0] == np.array([1,2])):
        roi[0] = 1064 - np.flip(smlROI[1])
    elif np.all(smlROI[0] == np.array([0,1])):
        pass
    roi[1] = smlROI[2]
    roi[0,0] = roi[0,0] - roiExp
    roi[0,1] = roi[0,1] + roiExp
    roi[1,0] = roi[1,0] - roiExp
    roi[1,1] = roi[1,1] + roiExp
    print(f'Cropping cube to {smlROI}, by small data notation, which is {roi} by our notation.')
    cropCubeFileByROI(exp,run,filename,roi)
    return
    # f = h5py.File(filename,'r')
    # scan_var = f['binVar_bins'][()]
    # bin_count = f['jungfrau1M_nEntries'][()]
    # imgs = f['jungfrau1M_data'][:,roi[0,0]:roi[0,1],roi[1,0]:roi[1,1]]
    # ROI = roi
    # i0 = f['ipm2__sum'][()]
    # i0_ipm3 = f['ipm3__sum'][()]
    # writeFile = h5py.File(filename[0:-3]+'_cropped.h5', 'w')
    # writeFile['scan_var'] = scan_var
    # writeFile['bin_count'] = bin_count
    # writeFile['imgs'] = imgs
    # writeFile['ROI'] = ROI
    # writeFile['i0'] = i0
    # writeFile['i0_ipm3'] = i0_ipm3
    # writeFile.close()
    # return
def cropRunCube(exp,run):
    #cubeDir = '/reg/data/drpsrcf/xpp/'+exp+'/scratch/hdf5/cube/'
    cubeDir = '/cds/data/drpsrcf/xpp/'+exp+'/scratch/hdf5/cube/'
    relevantFiles = []
    listOfFiles = os.listdir(cubeDir)
    for filename in listOfFiles:
        if f'{run:04}' in filename and 'cropped' not in filename:
            print(filename)
            relevantFiles.append(filename)
            if args.roi is not None:
                (x0, x1, y0, y1) = args.roi
                roi = np.zeros((2,2)).astype(int)
                roi[0,0] = x0
                roi[0,1] = x1
                roi[1,0] = y0
                roi[1,1] = y1
                cropCubeFileByROI(exp,run,filename,roi)
            else:
                cropCubeFile(exp,run,cubeDir+filename)
cropRunCube(exp,run)
