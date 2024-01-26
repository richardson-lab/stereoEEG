import sys
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import json
import tqdm

# Load data
root = sys.argv[1]
subject = sys.argv[2]

if len(sys.argv) < 3:
    print("Use: pacHist.py rootpath subject")
    sys.exit()

data = pd.read_csv(f'{root}/{subject}/{subject}_data.csv', header=None)
with open(f'{root}/{subject}/{subject}_param.json', 'r') as file:
    param = json.load(file)
    
if len(param['channel_id_labels']) == data.shape[1]:    
    print(f'Loaded {data.shape[0]} datapoints and {len(param.keys())} parameters for {data.shape[1]} channels')
else:
    print('Error loading data')
    
sample_rate = param['sample_rate']
tperm = param['tperm']
Nperm = param['Nperm']
Nch = data.shape[1]

sbin = round(param['tbin'] * sample_rate)
sstp = round(param['tstp'] * sample_rate)
rT = np.arange(0, data.shape[0]-sbin, sstp)
Nt = len(rT)
t = (rT + sbin / 2) / sample_rate
rP = np.column_stack((param['fP'][:-1], param['fP'][1:]))
rA = np.column_stack((param['fA'][:-1], param['fA'][1:]))
NP = rP.shape[0]
NA = rA.shape[0]
edges = np.arange(-np.pi, np.pi+param['delta'], param['delta'])
x = edges[:-1] + param['delta'] / 2
Nx = len(x)
PAC = np.zeros((Nch, NP, NA, Nt, Nx))
PACmi = np.zeros((Nch, NP, NA, Nt, 5))

for ich in tqdm(range(Nch), desc="Channel", leave=True):
    tmi = np.zeros((NP, NA, Nt))
    mrl = np.zeros((NP, NA, Nt))
    mu = np.zeros((NP, NA, Nt))
    tmip = np.zeros((NP, NA, Nt))
    mrlp = np.zeros((NP, NA, Nt))

    for iP in tqdm(range(NP), desc=f"Processing phase {iP}...", leave=False):
        bP, aP = butter(2, rP[iP, :] / (sample_rate / 2), btype='band')
        P = np.angle(hilbert(filtfilt(bP, aP, data[ich].to_numpy())))
        Pbin = np.zeros((sbin, Nt), dtype=int)
        for iT in range(Nt):
            t1 = rT[iT]
            t2 = rT[iT] + sbin
            cP = P[t1:t2]
            Pbin[:, iT] = np.digitize(cP, edges)
            
        for iA in tqdm(range(NA), desc=f"Processing amplitude {iA}...", leave=False):
            bA, aA = butter(2, rA[iA, :] / (sample_rate / 2), btype='band')
            A = np.abs(hilbert(filtfilt(bA, aA, data[ich].to_numpy())))
            for iT in range(Nt):
                t1 = rT[iT]
                t2 = rT[iT] + sbin
                cA = A[t1:t2]
                cAm = np.zeros(Nx)
                for jj in range(Nx):
                    cAm[jj] = np.mean(cA[Pbin[:, iT] == jj])
                PAC[ich, iP, iA, iT, :] = cAm
                cAm /= np.sum(cAm)
                cAm[cAm == 0] = 1e-10
                tmi[iP, iA, iT] = (np.log(Nx) + np.sum(cAm * np.log(cAm))) / np.log(Nx)
                r = np.sum(cA * np.exp(1j * cP))
                mrl[iP, iA, iT] = np.abs(r) / np.sum(cA)
                mu[iP, iA, iT] = np.angle(r)
            if Nperm > 1:
                MIperm = np.zeros((Nperm, Nt))
                Rperm = np.zeros((Nperm, Nt))
                for iperm in tqdm(range(Nperm), desc="Processing permutation test", leave=False):
                    tshift = tperm[0] + np.diff(tperm)[0] * np.random.rand()
                    nshift = round(tshift * sample_rate)
                    Ashift = np.roll(A, nshift)
                    for iT in range(Nt):
                        t1 = rT[iT]
                        t2 = rT[iT] + sbin
                        cA = Ashift[t1:t2]
                        cAm = np.zeros(Nx)
                        for jj in range(Nx):
                            cAm[jj] = np.mean(cA[Pbin[:, iT] == jj])
                        cAm /= np.sum(cAm)
                        cAm[cAm == 0] = 1e-10
                        MIperm[iperm, iT] = (np.log(Nx) + np.sum(cAm * np.log(cAm))) / np.log(Nx)
                        Rperm[iperm, iT] = np.abs(np.sum(cA * np.exp(1j * cP))) / np.sum(cA)

                n = np.sum(MIperm >= np.ones((Nperm, 1)) * np.squeeze(tmi[iP, iA, :]), axis=0)
                tmip[iP, iA, :] = (n + 1) / (Nperm + 1)
                n = np.sum(Rperm >= np.ones((Nperm, 1)) * np.squeeze(mrl[iP, iA, :]), axis=0)
                mrlp[iP, iA, :] = (n + 1) / (Nperm + 1)
                
    PACmi[ich] = np.concatenate([tmi[..., np.newaxis],
                             mrl[..., np.newaxis],
                             mu[..., np.newaxis],
                             tmip[..., np.newaxis],
                             mrlp[..., np.newaxis]], axis=3)
    
