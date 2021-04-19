import numpy as np
import pylab as pyl
import h5py

ldatfile = '/reg/data/drpsrcf/xpp/xpplv9818/scratch/hdf5/smalldata/xpplv9818_Run0020.h5'
f = h5py.File(ldatfile)
ttampl=f['tt']['AMPL'][()]
ids=ttampl>0.05
d=1e12*f['scan']['lxt_ttc'][()]
fltpos=f['tt']['FLTPOS'][()]
pyl.plot(fltpos[ids], d[ids], '.')
ttfit=np.polyfit(fltpos[ids], d[ids],1)
x = np.linspace(min(fltpos),max(fltpos),50)
#x = np.linspace(min(d),max(d),50)
pyl.plot(x, ttfit[0]*x+ttfit[1], "r")
pyl.show()
print ttfit

