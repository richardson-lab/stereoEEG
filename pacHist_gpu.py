import sys
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import json
from tqdm import tqdm
import h5py
import multiprocessing # remove?
from datetime import datetime
import cupy as cp
import cupyx.scipy.signal

def load_data_gpu(root, subject):
    data = pd.read_csv(f'{root}/{subject}/{subject}_data.csv', header=None)
    with open(f'{root}/{subject}/{subject}_param.json', 'r') as file:
        param = json.load(file)
    if len(param['channel_id_labels']) == data.shape[1]:    
        print(f'Loaded {data.shape[0]} datapoints and {len(param.keys())} parameters for {data.shape[1]} channels.')
        data_gpu = cp.asarray(data.values)
        return data_gpu, param
    else:
        print('Error loading data')
        return None, None
    
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
    rP = cp.asarray(np.column_stack((param['fP'][:-1], param['fP'][1:])))
    rA = cp.asarray(np.column_stack((param['fA'][:-1], param['fA'][1:])))
    NP = rP.shape[0]
    NA = rA.shape[0]
    edges = np.arange(-np.pi, np.pi+param['delta'], param['delta'])
    x = edges[:-1] + param['delta'] / 2
    Nx = len(x)
    return Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx

def perm_test_gpu(inputs):
    tperm, A, rT, sbin, sample_rate, Nt, Nx, Pbin, log_Nx = inputs
    
    # Convert inputs to GPU arrays
    A = cp.asarray(A)
    Pbin = cp.asarray(Pbin)
    MIperm = cp.zeros((1, Nt))
    
    tshift = tperm[0] + cp.diff(cp.asarray(tperm))[0] * cp.random.rand()
    nshift = round(tshift * sample_rate)
    Ashift = cp.roll(A, nshift)
    
    for iT in range(Nt):
        t1 = rT[iT]
        t2 = rT[iT] + sbin
        cA = Ashift[t1:t2]
        cAm = cp.array([cp.mean(cA[Pbin[:, iT] == jj]) for jj in range(1, Nx + 1)])
        cAm_sum = cp.sum(cAm)
        if cAm_sum > 0:
            cAm /= cAm_sum
            cAm[cAm == 0] = 1e-10
        MIperm[0, iT] = (log_Nx + cp.sum(cAm * cp.log(cAm))) / log_Nx
    
    # Convert MIperm back to a NumPy array if necessary for returning to CPU
    return MIperm.get()

def process_channel_gpu(data_gpu, param, shape):
    Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx = extract_params(param, shape)
    edges_gpu = cp.asarray(edges)
    
    tmi = cp.zeros((NP, NA, Nt))
    tmip = cp.zeros((NP, NA, Nt))
    log_Nx = cp.log(Nx)
    
    for iP in range(NP):
        bP, aP = cupyx.scipy.signal.butter(2, rP[iP, :] / (sample_rate / 2), btype='band', output='ba')
        P = cp.angle(cupyx.scipy.signal.hilbert(cupyx.scipy.signal.filtfilt(bP, aP, data_gpu)))
        Pbin_gpu = cp.transpose(cp.vstack([cp.digitize(P[t1:t1+sbin], edges_gpu) for t1 in rT]))
     
        for iA in range(NA):
            bA, aA = cupyx.scipy.signal.butter(2, rA[iA, :] / (sample_rate / 2), btype='band', output='ba')
            A = cp.abs(cupyx.scipy.signal.hilbert(cupyx.scipy.signal.filtfilt(bA, aA, data_gpu)))
            
            PAC = cp.zeros((Nch, NP, NA, Nt, Nx))
            for iT in range(Nt):
                cA = A[rT[iT]:rT[iT]+sbin]
                cAm = cp.array([cp.mean(cA[Pbin[:, iT] == jj+1]) for jj in range(Nx)])
                PAC[ich, iP, iA, iT, :] = cAm
                cAm /= cp.sum(cAm)
                cAm[cAm == 0] = 1e-10
                tmi[iP, iA, iT] = (log_Nx + cp.sum(cAm * cp.log(cAm))) / log_Nx
                    
            A_gpu = data_gpu[ich]  # Example: data for the current channel, already on GPU
            inputs = [(tperm, A_gpu, rT, sbin, sample_rate, Nt, Nx, Pbin_gpu, log_Nx) for _ in range(Nperm)]   

            # Prepare a batch of permutations to run in parallel
            # Example: Create an array of shape (Nperm, ...) where Nperm permutations
            # are to be processed in parallel
            permutations_gpu = cp.array([cp.roll(A_gpu, shift) for shift in shifts])

            # Call perm_test_gpu with batched inputs
            MIperm_gpu = perm_test_gpu(permutations_gpu, Pbin_gpu, other_params)

            
    PACmi_partial = cp.concatenate([tmi[..., cp.newaxis], tmip[..., cp.newaxis]], axis=3)
    
    # If you need to return the data to CPU memory:
    return PACmi_partial.get()  # This converts the result back to a NumPy array

def save_results_gpu(PACmi, root, subject):
    # Check if PACmi is a CuPy array and transfer it to CPU if necessary
    if isinstance(PACmi, cp.ndarray):
        PACmi_cpu = PACmi.get()  # Transfer from GPU to CPU
    else:
        PACmi_cpu = PACmi  # Assume it's already a CPU array if not a CuPy array
    
    # Proceed with the transposition and saving as before, using the CPU array
    PACmi_transposed = np.transpose(PACmi_cpu, (4, 3, 2, 1, 0))
    filename = f"{root}/{subject}/{subject}_PACmi.h5"
    with h5py.File(filename, 'w') as f:
        f.create_dataset('PACmi', data=PACmi_transposed)
        
def main_gpu():
    # Load data
    if len(sys.argv) < 3:
        print("Use: pacHist.py rootpath subject")
        sys.exit() 
    root = sys.argv[1]
    subject = sys.argv[2]
    
    # Loading data directly into GPU memory might not be straightforward for all formats
    # Convert the data to GPU array after loading if necessary
    data, param = load_data_gpu(root, subject)  # Assuming this function returns data ready for GPU
    
    # Extract parameters, assuming this function does not need GPU adaptation
    Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx = extract_params(param, data.shape)
    
    # Combine partial results, initiate PACmi on the GPU
    PACmi = cp.zeros((Nch, NP, NA, Nt, 2))
    
    with open("progress_log.txt", "w") as log_file:
        num_gpus = cupy.cuda.runtime.getDeviceCount()
        log_file.write(f"Number of GPUs available: {num_gpus}")
        log_file.flush() 
        
        for ich in range(Nch):
            log_file.write(f"Processing channel {ich+1}\n")
            log_file.flush()    
            start = datetime.now()
            
            channel_data_gpu = cp.asarray(data[ich]) if not isinstance(data[ich], cp.ndarray) else data[ich]
            PACmi[ich] = process_channel_gpu(channel_data_gpu, param, data.shape)
            
            end_time = datetime.now() - start
            log_file.write(f"Channel {ich+1} processing complete. Time elapsed: {end_time} \n \n")
            log_file.flush() 

    # Convert PACmi back to CPU for saving if necessary
    save_results_gpu(PACmi, root, subject)

if __name__ == "__main__":
    main_gpu()