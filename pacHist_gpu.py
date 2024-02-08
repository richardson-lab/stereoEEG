import sys
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import json
from tqdm import tqdm
import h5py
from datetime import datetime
import cupy as cp
import cupyx.scipy.signal


def load_data(root, subject):
    startLoad = datetime.now()
    data = pd.read_csv(f'{root}/{subject}/{subject}_data.csv', header=None)
    with open(f'{root}/{subject}/{subject}_param.json', 'r') as file:
        param = json.load(file)
    if len(param['channel_id_labels']) == data.shape[1]:
        print(
            f'Loaded {data.shape[0]} datapoints and {len(param.keys())} parameters for {data.shape[1]} channels. Time elapsed: {datetime.now() - startLoad}')
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
    rT = np.arange(0, datashape[0] - sbin, sstp)
    Nt = len(rT)
    t = (rT + sbin / 2) / sample_rate
    rP = np.column_stack((param['fP'][:-1], param['fP'][1:]))
    rA = np.column_stack((param['fA'][:-1], param['fA'][1:]))
    NP = rP.shape[0]
    NA = rA.shape[0]
    edges = np.arange(-np.pi, np.pi + param['delta'], param['delta'])
    x = edges[:-1] + param['delta'] / 2
    Nx = len(x)
    return Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx


def generate_ashifts(A, Nperm, tperm):
    interval_size = cp.diff(cp.array(tperm))[0]
    random_shifts = tperm[0] + interval_size * cp.random.rand(Nperm)
    N = A.size
    permuted_data = cp.empty((Nperm, N), dtype=A.dtype)
    for i, shift in enumerate(random_shifts):
        permuted_data[i, :] = cp.roll(A, int(shift))
    return permuted_data


def perm_test_gpu(Nperm, Nt, rT, sbin, Nx, log_Nx, a_shifts, Pbin_gpu):  # Updated this function
    MIperm = cp.zeros((Nperm, Nt))
    for iT in range(Nt):
        cA = a_shifts[:, rT[iT]:rT[iT] + sbin]
        cAm = cp.transpose(cp.array([cp.mean(cA[:, Pbin_gpu[:, iT] == jj], axis=1) for jj in range(1, Nx + 1)]))
        cAm_sum = cp.sum(cAm, axis=1)
        cAm = cAm / cAm_sum[:, cp.newaxis]
        MIperm[:, iT] = (log_Nx + cp.sum(cAm * cp.log(cAm), axis=1)) / log_Nx
    return MIperm


def process_channel_gpu(data_gpu, param, shape):
    Nch, sample_rate, tperm, Nperm, sbin, sstp, rT, Nt, t, rP, rA, NP, NA, edges, x, Nx = extract_params(param, shape)
    edges_gpu = cp.asarray(edges)
    tmi = cp.zeros((NP, NA, Nt))
    tmip = cp.zeros((NP, NA, Nt))
    log_Nx = cp.log(Nx)
    rT_gpu = cp.asarray(rT)

    for iP in tqdm(range(NP), desc="Channel progress", leave=False):
        bP, aP = butter(2, rP[iP, :] / (sample_rate / 2), btype='band')
        P = cp.angle(cupyx.scipy.signal.hilbert(cp.asarray(filtfilt(bP, aP, data_gpu.get()))))
        Pbin_gpu = cp.transpose(cp.vstack([cp.digitize(P[t1:t1 + sbin], edges_gpu) for t1 in rT]))

        for iA in tqdm(range(NA), desc=f"Phase {iP + 1} progress", leave=False):
            bA, aA = butter(2, rA[iA, :] / (sample_rate / 2), btype='band')
            A = cp.abs(cupyx.scipy.signal.hilbert(cp.asarray(filtfilt(bA, aA, data_gpu.get()))))

            for iT in range(Nt):
                cA = A[rT_gpu[iT]:rT_gpu[iT] + sbin]
                cAm = cp.array([cp.mean(cA[Pbin_gpu[:, iT] == jj + 1]) for jj in range(Nx)])
                cAm /= cp.sum(cAm)
                cAm[cAm == 0] = 1e-10
                tmi[iP, iA, iT] = (log_Nx + cp.sum(cAm * cp.log(cAm))) / log_Nx

            if Nperm > 0:
                a_shifts = generate_ashifts(A, Nperm, tperm)
                MIperm = perm_test_gpu(Nperm, Nt, rT, sbin, Nx, log_Nx, a_shifts, Pbin_gpu)
                n_gpu = cp.sum(MIperm >= cp.ones((Nperm, 1)) * cp.squeeze(tmi[iP, iA, :]), axis=0)
                tmip[iP, iA, :] = (n_gpu + 1) / (Nperm + 1)

    PACmi_partial = cp.concatenate([tmi[..., cp.newaxis], tmip[..., cp.newaxis]], axis=3)
    return PACmi_partial.get()


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
    if len(sys.argv) < 3:
        print("Use: pacHist.py rootpath subject")
        sys.exit()
    root = sys.argv[1]
    subject = sys.argv[2]
    data, param = load_data(root, subject)
    Nch, _, _, _, _, _, _, Nt, _, _, _, NP, NA, _, _, _ = extract_params(param, data.shape)
    PACmi = np.zeros((Nch, NP, NA, Nt, 2))
    with open("progress_log.txt", "w") as log_file:
        num_gpus = cp.cuda.runtime.getDeviceCount()
        log_file.write(f"Number of GPUs available: {num_gpus}\n")
        log_file.flush()
        for ich in range(Nch):
            log_file.write(f"Processing channel {ich + 1}\n")
            log_file.flush()
            start = datetime.now()
            channel_data_gpu = cp.asarray(data[ich]) if not isinstance(data[ich], cp.ndarray) else data[ich]
            PACmi[ich] = process_channel_gpu(channel_data_gpu, param, data.shape)
            end_time = datetime.now() - start
            log_file.write(f"Channel {ich + 1} processing complete. Time elapsed: {end_time} \n \n")
            log_file.flush()

    save_results_gpu(PACmi, root, subject)


if __name__ == "__main__":
    main_gpu()