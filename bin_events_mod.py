# import psana
import numpy as np

class binEvents(object):
  """Simple class to bin arrays based on a separate scalar value (usually delay).
  Assumes N bins between lb and ub, and stores N+2 elements in rows of _img.
  The averages for input x < min or x >= max are stored in 
  rows 0 and N+1, respectively.
  """
  def __init__(self, lb, ub, num_bins=10):
    self._min = lb
    self._max = ub
    self._img = None
    self._bin_width = float(ub-lb)/num_bins
    self._num_bins = num_bins + 2
    self._bin_count = np.zeros(self._num_bins)

  @classmethod
  def init_from_array(cls, a):
    delta = a[1]-a[0]
    amin = np.min(a)
    amax = np.max(a)
    mi = amin - delta/2
    ma = amax + delta/2
    n = len(a)
    return cls(mi, ma, n)

  def find_bin(self, x):
    """Finds corresponding bin for scalar value x.
    returns 0 or N+1 if the corresponding value is < _min or >= _max, respectively.
    """
    N = self._num_bins-2
    # real bins: 0 <= idx < N:
    idx = np.floor( (N)*(x-self._min)/(self._max - self._min) )  
    if idx < 0:
      return 0                # delays < min are stored in index 0
    if idx >= N:
      return N+1            # delays >= max are stored in index N+1
    return idx+1            # normal delays are in indices 1:N

  def update_bins(self, x, arr):
    bin_idx = self.find_bin(x)
    bin_idx = int(bin_idx)
    if arr is None:
      raise TypeError("'None' is not a valid input type.")
    else:
      if self._img is None:
        # initialize when we get the first element:
        if isinstance(arr, np.ndarray):
          arr_shape = [el for tuple in ((self._num_bins,), np.shape(arr)) for el in tuple]
          self._img = np.zeros(arr_shape, float)
        else:
          self._img = np.zeros((self._num_bins, np.size(arr) ), float)

      #n = self._bin_count[bin_idx]
      self._img[bin_idx,:] = self._img[bin_idx,:] + arr
      self._bin_count[bin_idx]+=1

  def bin_centers(self):
    num_bins = self._num_bins-2
    return np.linspace( self._min + self._bin_width/2, 
                        self._max - self._bin_width/2, 
                        num_bins)

  def bin_edges(self):
    num_bins = self._num_bins-2
    z = np.zeros(num_bins+1)
    for i in range(0,num_bins+1):
      z[i] = self._min + i*self._bin_width
    return z
    # np.linspace(self._min, self._max, self._num_bins-2)

