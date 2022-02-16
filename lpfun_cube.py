import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
from scipy.ndimage.filters import gaussian_filter
from mpi4py import MPI

comm = MPI.COMM_WORLD 
rank = comm.Get_rank()
size = comm.Get_size()

print(f'rank {rank} out of {size}')

def LPfun(data, varargin=None, justResults=False):
    if varargin is None:
        select_lambda_flag = True
    else:
        select_lambda_flag = False
        noo = varargin

    N = np.shape(data)[0]

    if N % 2 != 0:
        N -= 1

    x = data[:N, 1]
    t = data[:N, 0] - data[0, 0]

    deltat = t[1] - t[0]

    M = N//2

    # set up matrix from data (N-M)xM. It is backward prediction
    X = np.zeros((M,M), dtype='float')

    for i in range(M):
        X[:, i] = x[i+1:i+M+1]

    # computation of the (N-M)x(N-M) noonegativmatrix XX' and diagonalization

    XX = X@X.T
    d, U = np.linalg.eig(XX)
    # print('d', d)
    # d = d[::-1]
    d_length = len(d)

    dindex = np.argsort(d)
    d = d[dindex]
    U = U[:, dindex]
    no = np.arange(1, d_length+1) # change this to start from zero?

    if select_lambda_flag:
        plt.figure()
        plt.semilogy(no, d, 'o')
        plt.xlabel(r'$\lambda$')
        plt.ylabel('D')
        plt.show()
        noo = int(input('number of points lambda: '))

    # Filter eigenvalues
    E = np.zeros((M, M))

    for i in range(d_length-noo, d_length):
        E[i, i] = 1/d[i]

    xvector = x[:M]
    A = X.T @ U @ E @ U.T @ xvector
    # polynomial roots
    r = np.roots(np.concatenate((np.array([1]), -A.T))) # only l roots are significant being l rank of D matrix
    sindex = np.argsort(np.abs(r))
    s = r[sindex]
    ss = s.copy()
    ss[:] = s[::-1] # roots sorted in descending order

    b = np.log(np.abs(ss[:d_length])) / deltat
    w = np.angle(ss[:d_length]) / deltat

    P = np.sort(w)
    I = np.argsort(w)
    Z = b[I]

    # count number of w==0
    Nzeros = np.sum(w == 0)

    temp = (d_length-Nzeros)/2 + Nzeros
    temp - int(temp) + 1 # number between 1 and 2
    n_select = int(np.round(temp-int(temp)+1)) + int(temp) - 1
    # n_select = int(np.round((d_length-Nzeros)/2 + Nzeros))

    WW = np.abs(P[:n_select])
    BBB = Z[:n_select]

    # counting for positive damping constants

    Npos = np.sum(BBB >= 0)
    B1 = np.sort(BBB)
    J = np.argsort(BBB)
    W1 = WW[J]

    # sorted in descending order
    B2 = B1[::-1].copy()
    W2 = W1[::-1].copy()

    W = np.zeros(Npos)
    B = np.zeros(Npos)
    W[:Npos] = W2[:Npos]
    B[:Npos] = B2[:Npos]

    W_length = len(W)
    
    Xbar = np.zeros((N, W_length*2+1))
    i = np.arange(1, N+1)
    for j in range(W_length):
        Xbar[:, 2*j] = np.exp(-B[j]*(i-1)*deltat) * np.cos(W[j]*(i-1)*deltat)
        Xbar[:, 2*j+1] = -np.exp(-B[j]*(i-1)*deltat) * np.sin(W[j]*(i-1)*deltat)
    Xbar[:, W_length*2] = np.ones(N)

    AA = np.linalg.lstsq(Xbar, x, rcond=None)[0]

    C = np.zeros(W_length)
    fi = np.zeros(W_length)
    for i in range(W_length):
        if AA[2*i] == 0 and AA[2*i+1] == 0:
            C[i] = 0
            fi[i] = 0
        elif AA[2*i] == 0:
            fi[i] = np.sign(AA[2*i+1])*np.pi/2
            C[i] = np.abs(AA[2*i+1])
        elif AA[2*i+1] == 0:
            fi[i] = (1-np.sign(AA[2*i]))*np.pi/2
            C[i] = np.abs(AA[2*i])
        else:
            fi[i] = np.arctan2(AA[2*i+1], AA[2*i])
            C[i] = np.sqrt(AA[2*i+1]**2 + AA[2*i]**2)

    yy = np.zeros(N)
    zz = yy.copy()
    time = np.zeros((N, 3+W_length))
    time[:, 0] = t
    time[:, 1] = x
    
    timecc = np.zeros((N, W_length))
    Wmax = np.max(W)
    Wmax = Wmax / (2*np.pi)
    bi = B / (2*np.pi)
    ww = W / (2*np.pi)
    freq = np.arange(0, 1.5*Wmax, Wmax/1000)
    spect = np.zeros((len(freq), 1+W_length))
    spect[:, 0] = freq
    spect_1 = np.zeros((len(freq), W_length))

    totals = 0

    for i in range(W_length):
        timecc[:, i] = C[i] * np.exp(-B[i]*t)*np.cos(W[i]*t + fi[i])
        zz = zz + timecc[:, i]
        time[:, 3+i] = timecc[:, i]
        spect[:, 1+i] = C[i] * bi[i] / ((ww[i] - freq)**2 + bi[i]**2)
        totals = totals+spect[:, 1+i]
        spect_1[:, i] = spect[:, 1+i]

    time[:, 2] = zz
    Chi2 = np.sum((zz-x)**2)/N
    zz = zz + AA[2*W_length]

    # traces = h5py.File('LPtraces', 'w')
    # traces.create_dataset('t', data=t)
    # traces.create_dataset('y', data=timecc)
    # traces.create_dataset('each', data=spect_1)
    # traces.create_dataset('f', data=freq)
    # traces.create_dataset('sp', data=spect_1)

    if justResults==True:
        W = W / (2*np.pi)
        B = B / (2*np.pi)
        results = [W, B, fi, C]
        return results
    traces = {}
    traces['t'] = t
    traces['y'] = timecc
    traces['each'] = spect_1
    traces['f'] = freq
    traces['sp'] = spect_1
    traces['zz'] = zz
    traces['x'] = x
    traces['totals'] = totals

    W = W / (2*np.pi)
    B = B / (2*np.pi)
    results = [W, B, fi, C]

    return results, traces

def lpData(bins, normData,varargin=10):
    dataIn = np.vstack((bins,normData)).T
    #print(dataIn)
    try:
        results = LPfun(dataIn,varargin=varargin,justResults=True)
    except Exception as e:
        results = None
    return results

def lpDataArr(bins, normDataStack, varargin=10):
    imgShape = normDataStack.shape
    print(imgShape)
    resultStack = []
    for ii in range(0,imgShape[1]):
        results = lpData(bins,normDataStack[:,ii],varargin=varargin)
        resultStack.append(results)
        if ii%1000 == 0:
            print(f'Processed {ii} values.')
    return resultStack

def normDataByi0(imgs, i0):
    dataPoints = i0.shape[0]
    i0_flat = i0.flatten()
    normImg = imgs/i0[:,np.newaxis,np.newaxis]
    return normImg

def returnCubedata(filename):
    try:
        f = h5py.File(filename,'r')
    except Exception as e:
        print(f'Error opening cube file: {filename}, crop not saved for this cube.')
        return
    try:
        scan_var = f['binVar_bins'][()]
        imgs = f['jungfrau1M_data'][()]
        i0 = f['ipm2__sum'][()]
    except Exception as e:
        scan_var = f['scan_var']
        imgs = f['imgs']
        i0 = f['i0']
    return scan_var, i0, imgs



def gaussBlurData(data,sigma=7):
    for ii in range(0,data.shape[0]):
        data[ii] = gaussian_filter(data[ii], sigma=sigma)
    return data


if rank==0:

    def splitImgsAndSend(scan_var, normImg):
        dataShapes = []
        print(normImg.shape)
        splitNormImg = np.array_split(normImg, size, axis=1)
        for ii in range(1,size):
            data = splitNormImg[ii]
            dataShape = data.shape
            print(dataShape)
            dataShapes.append(dataShape)
            # comm.Send([np.array(dataShape), MPI.INT], dest=ii, tag=0)
            # comm.Send([data.flatten(), MPI.DOUBLE], dest=ii, tag=1)
            # comm.Send([scan_var, MPI.DOUBLE], dest=ii, tag=2)
            comm.send(np.array(dataShape), dest=ii, tag=0)
            comm.send(data.flatten(), dest=ii, tag=1)
            comm.send(scan_var, dest=ii, tag=2)
        return splitNormImg[0], splitNormImg[0].shape, dataShapes
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="Cube file name.", type=str)
    parser.add_argument("--saveFile", help="Save file name.", type=str,default=None)
    parser.add_argument('--zeroTime', help='At what index does the time signal start?', type=int,default=0)
    # parser.add_argument("--exp_name", help="specify the name of the experiment", type=str, default='xpply5120')

    args = parser.parse_args()
    filename = args.filename
    savefile = args.saveFile
    zeroTime = args.zeroTime
    if savefile is None:
        savefile = filename[0:-3]+'_lpfit.npz'
    # expname = args.exp_name

    scan_var, i0, imgs = returnCubedata(filename)
    imgs = imgs[zeroTime:]
    scan_var = scan_var[zeroTime:]
    i0 = i0[zeroTime:]
    normImg = normDataByi0(imgs, i0)
    #normImg = np.array(normImg,dtype=np.dtype(np.float32,align=True))
    imgs = None
    normImg = gaussBlurData(normImg,sigma=2)
    shapeOfimgs = normImg.shape
    normImg = normImg.reshape(shapeOfimgs[0],shapeOfimgs[1]*shapeOfimgs[2])
    data, dataShape, dataShapes = splitImgsAndSend(scan_var, normImg)
else:
    def receiveImg():
        # dataShape = np.empty(2,dtype=int)
        # comm.Recv(dataShape, source=0, tag=0)
        # data = np.empty(dataShape[0]*dataShape[1], dtype=np.float64)
        # comm.Recv(data, source=0, tag=1)
        # data = data.reshape(dataShape[0],dataShape[1])
        # scan_var = np.empty((dataShape[0],), dtype=np.float64)
        # comm.Recv(scan_var, source=0, tag=2)
        dataShape = comm.recv(source=0, tag=0)
        data = comm.recv(source=0, tag=1)
        data = data.reshape(dataShape[0],dataShape[1])
        scan_var = comm.recv(source=0, tag=2)
        return scan_var, data, dataShape
    scan_var, data, dataShape = receiveImg()

resultsStack = lpDataArr(scan_var, data, varargin=10)
#print(resultsStack)

if rank==0:
    def collectResults(resultsStack):
        receivedData = []
        receivedData.append(resultsStack)
        for ii in range(1,size):
            #curDataShape = dataShapes[ii-1]
            #recData = np.empty(curDataShape, dtype=np.float64)
            recData = comm.recv(source=ii, tag=4)
            receivedData.append(recData)
        return receivedData
    receivedData = collectResults(resultsStack)
    resultsStack = None
    receivedData = np.concatenate(receivedData,axis=0,dtype=object)
    receivedData = receivedData.reshape((shapeOfimgs[1],shapeOfimgs[2]))
    #writeFile = h5py.File(savefile, 'w')
    #dt = h5py.vlen_dtype(np.dtype('float32'))
    #dset = writeFile.create_dataset('lp_img', (shapeOfimgs[1],shapeOfimgs[2]), dtype=dt)
    #writeFile['zeroTime'] = zeroTime
    #writeFile['in_file'] = filename
    #dset[:,:] = receivedData
    np.savez(savefile, zeroTime=zeroTime, filename=filename,receivedData=receivedData)
else:
    def sendResults(resultsStack):
        comm.send(resultsStack, dest=0, tag=4)
        return
    sendResults(resultsStack)
