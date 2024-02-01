import sys
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import json
from tqdm import tqdm
import h5py
import multiprocessing
from datetime import datetime

## FUNCTIONS

def load_data(root, subject):
    startLoad = datetime.now()
    data = pd.read_csv(f'{root}/{subject}/{subject}_data.csv', header=None)
    with open(f'{root}/{subject}/{subject}_param.json', 'r') as file:
        param = json.load(file)
    if len(param['channel_id_labels']) == data.shape[1]:    
        print(f'Loaded {data.shape[0]} datapoints and {len(param.keys())} parameters for {data.shape[1]} channels. Time elapsed: {datetime.now() - startLoad}')
    else:
        print('Error loading data')
    return data, param

def extract_params(param, datashape):
    startParam = datetime.now()
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
    print(f"Loaded parameters. Time elapsed: {datetime.now() - startParam}")
    return Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx

def perm_test(inputs):
    startUnpack = datetime.now()
    tperm, A, rT, sbin, sample_rate, Nt, Nx, Pbin, log_Nx = inputs
    print(f"Unpacking inputs complete. Time elapsed: {datetime.now()-startUnpack}")
    MIperm = np.zeros((1,Nt));
    startrand = datetime.now()
    tshift = tperm[0] + np.diff(tperm)[0] * np.random.rand()
    print(f"Randomized tshift complete. Time elapsed: {datetime.now()-startrand}")
    startround = datetime.now()
    nshift = round(tshift * sample_rate)
    print(f"Rounded tshift complete. Time elapsed: {datetime.now()-startround}")
    startashift = datetime.now()
    Ashift = np.roll(A, nshift)
    print(f"Amplitude shift complete. Time elapsed: {datetime.now()-startashift}")
    startloop = datetime.now()
    for iT in range(Nt):
        t1 = rT[iT]
        t2 = rT[iT] + sbin
        cA = Ashift[t1:t2]
        cAm = np.array([np.mean(cA[Pbin[:, iT] == jj]) for jj in range(1, Nx + 1)])
        cAm_sum = np.sum(cAm)
        if cAm_sum > 0:
            cAm /= cAm_sum
            cAm[cAm == 0] = 1e-10
        MIperm[0,iT] = (log_Nx + np.sum(cAm * np.log(cAm))) / log_Nx
    print(f"Time loop complete. Time elapsed: {datetime.now()-startloop}")
    return MIperm

def process_channel(ich, data_numpy, param, shape):
    Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx = extract_params(param, shape)
    tmi = np.zeros((NP, NA, Nt))
    tmip = np.zeros((NP, NA, Nt))
    log_Nx = np.log(Nx)
    for iP in range(NP):
        startPhase = datetime.now()
        bP, aP = butter(2, rP[iP, :] / (sample_rate / 2), btype='band')
        P = np.angle(hilbert(filtfilt(bP, aP, data_numpy)))
        Pbin = np.transpose(np.vstack([np.digitize(P[t1:t1+sbin], edges) for t1 in rT]))
        print(f"Phase binning completed. Time elapsed: {datetime.now() - startPhase}")
     
        for iA in range(NA):
            startAmp = datetime.now()
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
            print(f"Amplitude binning completed. Time elapsed: {datetime.now() - startAmp}")
                  
            if Nperm > 1:
                Nperm = 5
                startInputs = datetime.now()
                inputs = [(tperm, A, rT, sbin, sample_rate, Nt, Nx, Pbin, log_Nx) for _ in range(Nperm)]
                print(f"Assigning permtest inputs completed. Time elapsed: {datetime.now() - startInputs}")
                startptest = datetime.now()
                MIperm = np.zeros((Nperm, Nt))
                for ipm in range(Nperm):
                      MIperm[ipm] = perm_test(inputs[ipm])
                print(f"Permutation test complete. Time elapsed: {datetime.now() - startInputs}")
                # with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                #     partial_perms = pool.map(perm_test, inputs)
                # MIperm = np.zeros((Nperm, Nt))
                # for iperm, partial_perm in enumerate(partial_perms):
                #     MIperm[iperm] = partial_perm
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
    
    # Combine partial results
    PACmi = np.zeros((Nch, NP, NA, Nt, 2))
    with open("progress_log.txt", "w") as log_file:
        log_file.write(f"Accessing {multiprocessing.cpu_count()} CPUs\n")
        log_file.flush() 
        for ich in range(Nch):
            log_file.write(f"Processing channel {ich+1}\n")
            log_file.flush()    
            start = datetime.now()
            PACmi[ich] = process_channel(ich, data[ich].to_numpy(), param, data.shape)
            log_file.write(f"Channel {ich+1} processing complete. Time elapsed: {datetime.now() - start} \n \n")
            log_file.flush() 

    # save_results(PACmi, root, subject)

if __name__ == "__main__":
    main()