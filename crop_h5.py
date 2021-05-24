# Img crop
# Copyright (C) 2021  Viktor Krapivin

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import h5py
import argparse
import numpy as np
import os
from scipy import ndimage

# parse arguments to script
parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="file to crop", type=str)
parser.add_argument("--outfile", help="filename to save crop as. Default is to append _crop to infile's name.", type=str)
parser.add_argument("--x", help="position of peak center, x coordinate. X and Y must both be specified, otherwise default method is used. Preempted by specifying ROI.", type=int)
parser.add_argument("--y", help="position of peak center, x coordinate", type=int)
parser.add_argument("--roi", help="coordinates of the edges of the roi rectangle in order of: x0, x1, y0, y1", type=int, nargs=4)
parser.add_argument("--crop_size", help="Side length of cropped square if one argument is given, or x side and then y side if two arguments are given. Default is 200 and this parameter is ignored if roi is specified.", type=int, nargs='+')
parser.add_argument("--method", help='How to identify peak location. The current methods supported are min, max, gmin, gmax. The last two parameters correspond to the Gaussian fourier filtered min and max. Methods are applied to average image difference between on-off during the last 20 scans. To modify the images this is applied to use the next argument. Default is gmin.', type=str)
parser.add_argument("--method_images", help='How to choose image on which to peak find. Options are: default last 20 images of (on-off), late_early (last 20 images in scan minus first 20 images in scan, both laser on). To specify different number of images specify the method and append an image count with an underscope. Ie. if you would like to use last 30 images uses default_30 or late_early_30.', type=str)
parser.add_argument("--num_avg_img", help="Number of images to average.", type=int,default=20)
args = parser.parse_args()



numofOfImagesToAverage = args.num_avg_img
# print(args)
# print(dir(args))
# for ii in dir(args):
#     exec('print(args.'+ii+')')

usingCoords = False
if (args.roi is not None):
    print('cropping with roi specified')
elif (args.x is None) or (args.y is None):
    print('cropping around peak found with automatic peak finding')
else:
    usingCoords = True
    print('cropping using provided coordinates.')

h5filename = args.infile
print('opening file ' + h5filename)
cubed = h5py.File(h5filename,'r')
i0_laser_on = cubed['ipm2_sum'][:].flatten()
i0_laser_off = cubed['ipm2_sum_off'][:].flatten()
i0_laser_on_extra = cubed['ipm3_sum'][:].flatten()
i0_laser_off_extra = cubed['ipm3_sum_off'][:].flatten()
detector_laser_on = cubed['cspad'][:]
detector_laser_off = cubed['cspad_off'][:]
scan_axis =  cubed['bins'][:].flatten()
bin_counts = cubed['nEntries'][:].flatten()
bin_counts_off = cubed['nEntries_off'][:].flatten()
cubed.close()
print('read cube, processing now')

h5fileout = args.outfile
if h5fileout is None:
    h5fileout = h5filename[0:len(h5filename)-3] + "_cropped.h5"
print('creating h5 crop file: ' + h5fileout)
h5f_out = h5py.File(h5fileout, 'w')

h5f_out.create_dataset('I0_laser_on_integrated', data=i0_laser_on)
h5f_out.create_dataset('I0_laser_on_integrated_off', data=i0_laser_off)
h5f_out.create_dataset('bin_counts', data=bin_counts)
h5f_out.create_dataset('bin_counts_off', data=bin_counts_off)
h5f_out.create_dataset('scan_var_position', data=scan_axis)


# impose image boundaries on ROI
def impose_image_boundaries_on_roi(roiArr, imageshape):
    roiCoordsTemp = np.copy(roiArr)
    # make sure coordinates are not greater than image shape
    for ii in range(0,roiCoordsTemp.shape[0]):
        if roiCoordsTemp[ii] < 0:
            roiCoordsTemp[ii] = 0
    if roiCoordsTemp[0] > imageshape[0]:
        roiCoordsTemp[0] = imageshape[0]
    if roiCoordsTemp[1] > imageshape[0]:
        roiCoordsTemp[1] = imageshape[0]
    if roiCoordsTemp[2] > imageshape[1]:
        roiCoordsTemp[2] = imageshape[1]
    if roiCoordsTemp[3] > imageshape[1]:
        roiCoordsTemp[3] = imageshape[1]
    return roiCoordsTemp
# take the crop size arg and process to get crop x and crop y
def crop_size_to_ROI(cropsize_arg):
    crop_size_x = 200
    crop_size_y = 200
    if cropsize_arg is None:
        pass
    elif len(cropsize_arg) == 1:
        crop_size_x = cropsize_arg[0]
        crop_size_y = cropsize_arg[0]
    elif len(cropsize_arg) == 2:
        crop_size_x = cropsize_arg[0]
        crop_size_y = cropsize_arg[1]
    return (crop_size_x, crop_size_y)
# use crop size and impose image boundaries to generate roi tuple.
def roi_coords_to_crop_coords(peakX, peakY, crop_size_arg, detector_image_shape):
    (crop_size_x, crop_size_y) = crop_size_to_ROI(crop_size_arg)
    roiCoords = np.array([peakX - np.round(crop_size_x/2.0), peakX + np.round(crop_size_x/2.0), peakY - np.round(crop_size_y/2.0), peakY + np.round(crop_size_y/2.0)])
    roiCoords = impose_image_boundaries_on_roi(roiCoords, detector_image_shape)
    return roiCoords
def save_det_images(h5f_out, peakX, peakY, crop_size_arg, detector_laser_on, detector_laser_off):
    roiCoords = roi_coords_to_crop_coords(peakX, peakY, crop_size_arg, detector_laser_on.shape[1:3])
    (x0, x1, y0, y1) = tuple(roiCoords.astype(int))
    h5f_out.create_dataset('detector_laser_on', data=detector_laser_on[:,x0:x1,y0:y1])
    h5f_out.create_dataset('detector_laser_off', data=detector_laser_off[:,x0:x1,y0:y1])
    h5f_out.create_dataset('roi', data=roiCoords)
    return (x0, x1, y0, y1)


if args.roi is not None:
    (x0, x1, y0, y1) = args.roi
    h5f_out.create_dataset('detector_laser_on', data=detector_laser_on[:,x0:x1,y0:y1])
    h5f_out.create_dataset('detector_laser_off', data=detector_laser_off[:,x0:x1,y0:y1])
    print(f'done cropping. Used user-specified ROI: ({x0}, {x1}, {y0}, {y1})')
    roiCoords = np.array([x0, x1, y0, y1])
    h5f_out.create_dataset('roi', data=roiCoords)
elif usingCoords == True:
    peakX = args.x
    peakY = args.y
    (x0, x1, y0, y1) = save_det_images(h5f_out, peakX, peakY, args.crop_size, detector_laser_on, detector_laser_off)
    print(f'done cropping. Generated ROI: ({x0}, {x1}, {y0}, {y1}) from user specified peak center and crop size')
else:
    on_img = detector_laser_on/i0_laser_on[:,None,None]
    off_img = detector_laser_off/i0_laser_off[:,None,None]
    # find bin_counts that are above the mean
    bin_count_amean = np.where(bin_counts >= np.mean(bin_counts))[0]
    bin_count_end_ind = int(bin_count_amean[-20].astype(int))
    if args.method_images is not None:
        numofOfImagesToAverage = int(args.method_images.split('_')[-1:][0])
    bin_count_start_ind = bin_count_end_ind - numofOfImagesToAverage

    process_img_to_use = None
    if args.method_images is not None and 'late_early' in args.method_images:
        process_img_to_use = np.sum(on_img[bin_count_start_ind:bin_count_end_ind,:,:] - off_img[bin_count_start_ind:bin_count_end_ind,:,:],axis=0)
    else:
        process_img_to_use = np.sum(on_img[bin_count_start_ind:bin_count_end_ind,:,:] - off_img[bin_count_start_ind:bin_count_end_ind,:,:],axis=0)
    if args.method == 'min':
        pass
    elif args.method == 'max':
        pass
    else: # default gmin, gmax included as well.
        process_img_to_use = ndimage.gaussian_filter(np.sum(on_img[bin_count_start_ind:bin_count_end_ind,:,:] - off_img[bin_count_start_ind:bin_count_end_ind,:,:],axis=0), 10)
    
    
    if args.method is not None and 'max' in args.method:
        peakX, peakY = np.unravel_index(np.nanargmax(process_img_to_use), process_img_to_use.shape)
    else: #includes 'min' 
        print(np.unravel_index(np.nanargmin(process_img_to_use), process_img_to_use.shape))
        peakX, peakY = np.unravel_index(np.nanargmin(process_img_to_use), process_img_to_use.shape)
    (x0, x1, y0, y1) = save_det_images(h5f_out, peakX, peakY, args.crop_size, detector_laser_on, detector_laser_off)
    print(f'done cropping. Automatically generated ROI as: ({x0}, {x1}, {y0}, {y1}) from supplied algorithm arguments.')
    




h5f_out.close()