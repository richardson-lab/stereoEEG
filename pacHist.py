import sys
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import json
from tqdm import tqdm
import h5py
import multiprocessing
from functools import partial
import cProfile

## FUNCTIONS

def load_data(root, subject):
    data = pd.read_csv(f'{root}/{subject}/{subject}_data.csv', header=None)
    with open(f'{root}/{subject}/{subject}_param.json', 'r') as file:
        param = json.load(file)
    if len(param['channel_id_labels']) == data.shape[1]:    
        print(f'Loaded {data.shape[0]} datapoints and {len(param.keys())} parameters for {data.shape[1]} channels')
    else:
        print('Error loading data')
    return data, param

def extract_params(param, datashape):
    Nch = datashape[1]
    sample_rate = param['sample_rate']
    tperm = param['tperm']
    Nperm = param['Nperm']
    sbin = round(param['tbin'] * sample_rate)
    sstp = round(param['tstp'] * sample_rate)
    rT = np.arange(0, datashape[0]-sbin, sstp)
    Nt = len(rT)
    t = (rT + sbin / 2) / sample_rate
    rP = np.column_stack((param['fP'][:-1], param['fP'][1:]))
    rA = np.column_stack((param['fA'][:-1], param['fA'][1:]))
    NP = rP.shape[0]
    NA = rA.shape[0]
    edges = np.arange(-np.pi, np.pi+param['delta'], param['delta'])
    x = edges[:-1] + param['delta'] / 2
    Nx = len(x)
    return Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx

def process_channel(ich, data, param):
    Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx = extract_params(param, data.shape)

    tmi = np.zeros((NP, NA, Nt))
    tmip = np.zeros((NP, NA, Nt))

    data_numpy = data[ich].to_numpy()
    log_Nx = np.log(Nx)

    for iP in tqdm(range(NP-15), desc=f'Channel {ich} Phase Progress'):
        bP, aP = butter(2, rP[iP, :] / (sample_rate / 2), btype='band')
        P = np.angle(hilbert(filtfilt(bP, aP, data_numpy)))
        Pbin = np.transpose(np.vstack([np.digitize(P[t1:t1+sbin], edges) for t1 in rT]))
            
        for iA in range(NA-15):
            bA, aA = butter(2, rA[iA, :] / (sample_rate / 2), btype='band')
            A = np.abs(hilbert(filtfilt(bA, aA, data_numpy)))
            
            PAC = np.zeros((Nch, NP, NA, Nt, Nx))
            for iT in range(Nt):
                cA = A[rT[iT]:rT[iT]+sbin]
                cAm = np.array([np.mean(cA[Pbin[:, iT] == jj+1]) for jj in range(Nx)])
                PAC[ich, iP, iA, iT, :] = cAm
                cAm /= np.sum(cAm)
                cAm[cAm == 0] = 1e-10
                tmi[iP, iA, iT] = (log_Nx + np.sum(cAm * np.log(cAm))) / log_Nx
                
            if Nperm > 1:
                MIperm = np.zeros((Nperm, Nt))
                for iperm in range(Nperm-195):
                    tshift = tperm[0] + np.diff(tperm)[0] * np.random.rand()
                    nshift = round(tshift * sample_rate)
                    Ashift = np.roll(A, nshift)
                    
                    for iT in range(Nt):
                        t1 = rT[iT]
                        t2 = rT[iT] + sbin
                        cA = Ashift[t1:t2]
                        cAm = np.zeros(Nx)
                        for jj in range(Nx):
                            cAm[jj] = np.mean(cA[Pbin[:, iT] == jj+1])
                        cAm /= np.sum(cAm)
                        cAm[cAm == 0] = 1e-10
                        MIperm[iperm, iT] = (log_Nx + np.sum(cAm * np.log(cAm))) / log_Nx
                n = np.sum(MIperm >= np.ones((Nperm, 1)) * np.squeeze(tmi[iP, iA, :]), axis=0)
                tmip[iP, iA, :] = (n + 1) / (Nperm + 1)
    PACmi_partial = np.concatenate([tmi[..., np.newaxis],tmip[..., np.newaxis]],axis=3)
    return PACmi_partial            
    
def save_results(PACmi, root, subject):
    PACmi_transposed = np.transpose(PACmi, (4, 3, 2, 1, 0))
    filename = f"{root}/{subject}/{subject}_PACmi.h5"
    with h5py.File(filename, 'w') as f:
        f.create_dataset('PACmi', data=PACmi_transposed)

## MAIN FUNCTION

def main():

    # Load data
    if len(sys.argv) < 3:
        print("Use: pacHist.py rootpath subject")
        sys.exit() 
    root = sys.argv[1]
    subject = sys.argv[2]
    data, param = load_data(root, subject)
    Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx = extract_params(param, data.shape)
    
    # Set up and execute multiprocessing
    process_with_args = partial(process_channel, data=data, param=param)
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    partial_results = list(tqdm(pool.imap(process_with_args, range(Nch)), total=Nch, leave=True)) # local
    # partial_results = pool.map(process_with_args, range(Nch)) # HPC
    pool.close()
    pool.join()

    # Combine partial results
    PACmi = np.zeros((Nch, NP, NA, Nt, 2))
    for ich, partial_result in enumerate(partial_results):
        PACmi[ich] = partial_result

    save_results(PACmi, root, subject)

if __name__ == "__main__":
    cProfile.run(main())