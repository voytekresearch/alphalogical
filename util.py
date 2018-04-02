import os
import csv
import h5py
import warnings

import numpy as np

from copy import deepcopy
from collections import defaultdict

from scipy.io import loadmat
from scipy.signal import butter, lfilter
from scipy.signal import hilbert

from fakespikes.util import bin_times
from fakespikes.util import create_times

# --------------------------------------------------------
# Setup cache
from tempfile import mkdtemp
from joblib import Memory
if "ALPHA_CACHEDIR" in os.environ:
    cachedir = os.environ["ALPHA_CACHEDIR"]
else:
    cachedir = mkdtemp()
    print("ALPHA_CACHEDIR was not preset.")
    print("ALPHA_CACHEDIR set to {}".format(cachedir))

memory = Memory(cachedir=cachedir, verbose=0)

# --------------------------------------------------------


def find_closest_idx(t, times):
    times = np.asarray(times)
    idx = (np.abs(times - t)).argmin()

    return idx


def find_closest(t, times):
    idx = find_closest_idx(t, times)
    return times[idx]


def load_good_smith_files(data_path):
    lfp_files = [
        "Bo130408_s6ae_fixblank_active_0001_converted_ns2.mat",
        "Bo130408_s6ae_fixblank_active_0002_converted_ns2.mat",
        "Bo130409_s7ae_fixblank_active_0003_converted_ns2.mat",
        "Wi130116_s51ae_fixblank_active_0001_converted_ns2.mat",
        "Bo130404_s4ae_fixblank_active_0002_converted_ns2.mat",
        "Bo130405_s5ae_fixblank_active_0001_converted_ns2.mat",
        "Bo130405_s5ae_fixblank_active_0002_converted_ns2.mat",
        "Bo130405_s5ae_fixblank_active_0003_converted_ns2.mat",
        "Bo130418_s12ae_fixblank_active_0002_converted_ns2.mat",
        "Wi121219_s43ae_fixblank_active_0001_converted_ns2.mat",
        "Wi121219_s43ae_fixblank_active_0002_converted_ns2.mat",
        "Wi130129_s55ae_fixblank_active_0001_converted_ns2.mat",
        "Wi130129_s55ae_fixblank_active_0002_converted_ns2.mat",
        "Wi130205_s58ae_fixblank_active_0001_converted_ns2.mat",
        "Wi130205_s58ae_fixblank_active_0002_converted_ns2.mat",
        "Wi130207_s59ae_fixblank_active_0001_converted_ns2.mat",
        "Wi130207_s59ae_fixblank_active_0002_converted_ns2.mat",
        "Wi130207_s59ae_fixblank_active_0003_converted_ns2.mat",
        "Wi130208_s60ae_fixblank_active_0001_converted_ns2.mat",
        "Wi130211_s61ae_fixblank_active_0001_converted_ns2.mat",
        "Wi130212_s62ae_fixblank_active_0001_converted_ns2.mat"
    ]

    # Generate rate file names based on the good list
    rate_files = []
    for fi in lfp_files:
        fi_name = os.path.splitext(fi)[0]
        fi_name = fi_name.split("_")[:5]
        fi_name = "_".join(fi_name)

        new_name = "{}.mat".format(fi_name)

        rate_files.append(new_name)

    # Generate seg file names based on the good list
    seg_files = []
    for fi in lfp_files:
        fi_name = os.path.splitext(fi)[0]
        new_name = "{}_segments.csv".format(fi_name)

        seg_files.append(new_name)

    # Add pathing
    lfp_files = [os.path.join(data_path, fi) for fi in lfp_files]
    rate_files = [os.path.join(data_path, fi) for fi in rate_files]
    seg_files = [os.path.join(data_path, fi) for fi in seg_files]

    return lfp_files, rate_files, seg_files


def load_smith_foof_results(name):
    results = []
    with open(name, 'r') as fi:
        reader = csv.reader(fi, delimiter=",")
        for i, row in enumerate(reader):
            if i == 0:
                header = row
                continue

            # Convert elements to the right type
            # ['m', 'c', 'i', 'j', 'center', 'power', 'bw']
            # [int, int,  int, int, float, ...]
            m, c, i, j, center, power, bw = row
            m = int(m)
            c = int(c)
            i = int(i)
            j = int(j)
            center = float(center)
            power = float(power)
            bw = float(bw)

            typed_row = (m, c, i, j, center, power, bw)
            results.append(typed_row)

    return header, results


@memory.cache
def load_smith_rates(name, channel_num, t_range, dt):
    """Load Smith rate data."""

    # Load .mat from Smith
    data = loadmat(name)['nevDat']

    # Parse it
    event_times = data[:, 2]
    rate_window = 1e-3

    # Get channel data
    channel_mask = channel_num == data[:, 0]
    event_times = data[channel_mask, 2]

    # Create rates
    times, rates = bin_times(event_times, t_range, rate_window)

    return times, rates


@memory.cache
def load_smith_lfps(name, channel_num):
    """Load Smith LFP data.
    
    NOTE: .ns2 files must first be converted to .mat
    """

    # Load data
    fi = h5py.File(name)

    # Parse time
    fs = float(fi['Fs'].value)
    n_samples = float(fi['nSamples'].value)
    T = n_samples / fs

    times = create_times(T, 1 / fs)

    # Load data
    data = fi['data'].value

    return fs, times, data[:, channel_num]


# --
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')

    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)

    return y


@memory.cache
def hilbert_power(x):
    return np.abs(hilbert(x))


def write_analysis(name, results, mode="w"):
    name = "{}.csv".format(name)

    # If the file exists, don't append another
    # header....
    header = True
    if os.path.isfile(name) and mode == "a":
        header = False

    # Write that data!
    keys = sorted(results.keys())
    with open(name, mode=mode) as fi:
        writer = csv.writer(fi, delimiter=",")

        if header:
            writer.writerow(keys)

        writer.writerows(zip(* [results[key] for key in keys]))


def shift_bursts(bursts, shift):
    shifted = []
    for b in bursts:
        shifted.append([i + shift for i in b])

    return shifted


def shift_extrema(extrema, shift):
    extrema = [ex + shift for ex in extrema]
    return extrema


def analyze_bursts(analysis_code, bursts, data, fs=1000):
    if len(bursts) == 0:
        return None

    results = defaultdict(list)
    for k, b in enumerate(bursts):
        b = np.asarray(b)

        # -
        # Stats
        n = k
        length = float(len(b))
        time = length / fs
        mean = np.mean(data[b])
        med = np.median(data[b])
        sd = np.std(data[b])

        # -
        results["analysis_code"].append(analysis_code)
        results["count"].append(k)
        results["length"].append(length)
        results["time"].append(time)

        results["mean"].append(mean)
        results["med"].append(med)
        results["sd"].append(sd)

    return results


def analyze_extrema(analysis_code, extrema, data, n_sample=10, fs=1000):
    if len(extrema) == 0:
        return None

    results = defaultdict(list)
    for k, ex in enumerate(extrema):
        b = range(ex - int(n_sample / 2), ex + int(n_sample / 2))
        b = np.asarray(b)

        # -
        # Stats
        n = k
        mean = np.mean(data[b])
        med = np.median(data[b])
        sd = np.std(data[b])

        # -
        results["analysis_code"].append(analysis_code)
        results["count"].append(k)

        results["mean"].append(mean)
        results["med"].append(med)
        results["sd"].append(sd)

    return results


def find_bursts(x, start=4, end=2):
    x = deepcopy(x)
    #     x -= np.median(x)
    #     x /= np.std(x)

    bursts = []
    burst = []
    b = False
    for i, xt in enumerate(x):
        # start of a burst
        if xt >= start:
            b = True

        if b:
            # Still bursting?
            if xt >= end:
                burst.append(i)
            else:
                # we're below the thresh. store
                bursts.append(deepcopy(burst))

                # and then reset
                burst = []
                b = False

    return bursts


def find_notbursts(x, zero_thresh):
    n = x.shape[0]

    # Find zeros
    candidates = np.zeros_like(x, dtype=np.bool)
    for i, xt in enumerate(x):
        if xt < zero_thresh:
            candidates[i] = True

    # Convert zeros into notburst indices
    notbursts = []
    notburst = []
    b = False
    for i in range(n):
        # notburst at i?
        c = candidates[i]

        # Still notbursting?
        if b and c:
            notburst.append(i)
        # notburst onset?
        elif c:
            notburst = [
                i,
            ]
            b = True
        # notburst over.
        elif b:
            # Store that notbursts index
            notbursts.append(deepcopy(notburst))

            # and reset
            notburst = []
            b = False
        else:
            notburst = []
            b = False

    return notbursts


def find_bursts2(x, burst_thresh, zero_thresh, history=200):
    # Find zeros
    zeros = np.zeros_like(x, dtype=np.bool)
    for i, xt in enumerate(x):
        if xt < zero_thresh:
            zeros[i] = True

    # Find candidate busts
    candidates = np.zeros_like(x, dtype=np.bool)
    for i, xt in enumerate(x):
        if xt >= burst_thresh:
            candidates[i] = True

    # Define bursts, checking that zeros preceed bursts
    n = x.shape[0]

    bursts = []
    burst = []
    b = False
    for i in range(history, n):
        # History is filled with zeros and not bursts?
        z_last = zeros[history - i:i].sum()
        c_last = candidates[history - i:i].sum()
        if (z_last > 0) and np.allclose(c_last, 0.0):
            z_last = True
        else:
            z_last = False

        # burst at i?
        c = candidates[i]

        # Still bursting?
        if b and c:
            burst.append(i)
        # Burst onset?
        elif c and z_last:
            burst = [
                i,
            ]
            b = True
        # Burst over.
        elif b:
            # Store that bursts index
            bursts.append(deepcopy(burst))

            # and reset
            burst = []
            b = False
        else:
            burst = []
            b = False

    return bursts


def find_extrema(x_filt, Fs, f_range, boundary=None, first_extrema='peak'):
    """
    Identify peaks and troughs in a time series.

    NOTE: This code is based on: https://github.com/voytekresearch/neurodsp

    That code was ported, sans filtering the component, to py2 to match
    the rest of the alphalogical code. 
    
    Parameters
    ----------
    x_filt : array-like 1d
        A filtered voltage time series
    Fs : float
        sampling rate
    f_range : (low, high), Hz
        frequency range for narrowband signal of interest,
        used to find zerocrossings of the oscillation
    boundary : int
        number of samples from edge of recording to ignore
    first_extrema: str or None
        if 'peak', then force the output to begin with a peak and end in a trough
        if 'trough', then force the output to begin with a trough and end in peak
        if None, force nothing
    filter_fn : filter function, `filterfn(x, Fs, pass_type, f_lo, f_hi, remove_edge_artifacts=True)
        Must have the same API as neurodsp.filter
    filter_kwargs : dict
        keyword arguments to the filter_fn

    Returns
    -------
    Ps : array-like 1d
        indices at which oscillatory peaks occur in the input signal x
    Ts : array-like 1d
        indices at which oscillatory troughs occur in the input signal x

    Notes
    -----
    This function assures that there are the same number of peaks and troughs
    if the first extrema is forced to be either peak or trough.
    """

    # Default boundary value as 1 cycle length
    if boundary is None:
        boundary = int(np.ceil(Fs / float(f_range[0])))

    # Find rising and falling zerocrossings
    zeroriseN = _fzerorise(x_filt)
    zerofallN = _fzerofall(x_filt)

    # Compute number of peaks and troughs
    if zeroriseN[-1] > zerofallN[-1]:
        P = len(zeroriseN) - 1
        T = len(zerofallN)
    else:
        P = len(zeroriseN)
        T = len(zerofallN) - 1

    # Calculate peak samples
    Ps = np.zeros(P, dtype=int)
    for p in range(P):
        # Calculate the sample range between the most recent zero rise
        # and the next zero fall
        mrzerorise = zeroriseN[p]
        nfzerofall = zerofallN[zerofallN > mrzerorise][0]
        # Identify time fo peak
        Ps[p] = np.argmax(x_filt[mrzerorise:nfzerofall]) + mrzerorise

    # Calculate trough samples
    Ts = np.zeros(T, dtype=int)
    for tr in range(T):
        # Calculate the sample range between the most recent zero fall
        # and the next zero rise
        mrzerofall = zerofallN[tr]
        nfzerorise = zeroriseN[zeroriseN > mrzerofall][0]
        # Identify time of trough
        Ts[tr] = np.argmin(x_filt[mrzerofall:nfzerorise]) + mrzerofall

    # Remove peaks and troughs within the boundary limit
    Ps = Ps[np.logical_and(Ps > boundary, Ps < len(x_filt) - boundary)]
    Ts = Ts[np.logical_and(Ts > boundary, Ts < len(x_filt) - boundary)]

    # Force the first extrema to be as desired
    # Assure equal # of peaks and troughs
    if first_extrema == 'peak':
        if Ps[0] > Ts[0]:
            Ts = Ts[1:]
        if Ps[-1] > Ts[-1]:
            Ps = Ps[:-1]
    elif first_extrema == 'trough':
        if Ts[0] > Ps[0]:
            Ps = Ps[1:]
        if Ts[-1] > Ps[-1]:
            Ts = Ts[:-1]
    elif first_extrema is None:
        pass
    else:
        raise ValueError('Parameter forcestart is invalid')

    return Ps, Ts


def _fzerofall(data):
    """Find zerocrossings on falling edge of a filtered signal"""
    pos = data > 0
    return (pos[:-1] & ~pos[1:]).nonzero()[0]


def _fzerorise(data):
    """Find zerocrossings on rising edge of a filtered signal"""
    pos = data < 0
    return (pos[:-1] & ~pos[1:]).nonzero()[0]


def find_zerox(x, Ps, Ts):
    """
    Find zerocrossings within each cycle after peaks and troughs are identified.
    A rising zerocrossing occurs when the voltage crosses
    midway between the trough voltage and subsequent peak voltage.
    A decay zerocrossing is defined similarly.
    If this voltage is crossed at multiple times, the temporal median is taken
    as the zerocrossing.

    NOTE: This code is based on: https://github.com/voytekresearch/neurodsp

    Parameters
    ----------
    x : array-like 1d
        voltage time series
    Ps : numpy arrays 1d
        time points of oscillatory peaks
    Ts : numpy arrays 1d
        time points of osillatory troughs

    Returns
    -------
    zeroxR : array-like 1d
        indices at which oscillatory rising zerocrossings occur
    zeroxD : array-like 1d
        indices at which oscillatory decaying zerocrossings occur

    Notes
    -----
    * Sometimes, due to noise in estimating peaks and troughs when the oscillation
    is absent, the estimated peak might be lower than an adjacent trough. If this
    occurs, the rise and decay zerocrossings will be set to be halfway between
    the peak and trough.
    """

    # Calculate the number of rises and decays
    if Ps[0] < Ts[0]:
        N_rises = len(Ps) - 1
        N_decays = len(Ts)
        idx_bias = 0
    else:
        N_rises = len(Ps)
        N_decays = len(Ts) - 1
        idx_bias = 1

    # Find zerocrossings for rise
    zeroxR = np.zeros(N_rises, dtype=int)
    for i in range(N_rises):
        x_temp = np.copy(x[Ts[i]:Ps[i + 1 - idx_bias] + 1])
        x_temp -= (x_temp[0] + x_temp[-1]) / 2.

        # Catch if rise is actually a net decay
        if x_temp[0] > x_temp[-1]:
            zeroxR[i] = Ts[i] + int(len(x_temp) / 2.)
        else:
            try:
                zeroxR[i] = Ts[i] + int(np.median(_fzerorise(x_temp)))
            except:
                warnings.warn(
                    'Error when estimating rising zerocrossing after trough ' +
                    str(i) + ' at sample ' + str(Ts[i]) +
                    '. Therefore, the zerocrossing has been set to halfway between the two extrema.'
                    +
                    ' This is potentially due to the two extrema being too close together'
                )
                zeroxR[i] = Ts[i] + int(len(x_temp) / 2.)

    # Find zerocrossings for decays
    zeroxD = np.zeros(N_decays, dtype=int)
    for i in range(N_decays):
        x_temp = np.copy(x[Ps[i]:Ts[i + idx_bias] + 1])
        x_temp -= (x_temp[0] + x_temp[-1]) / 2.

        # Catch if the decay period is actually a net rise
        if x_temp[0] < x_temp[-1]:
            zeroxD[i] = Ps[i] + int(len(x_temp) / 2.)
        else:
            try:
                zeroxD[i] = Ps[i] + int(np.median(_fzerofall(x_temp)))
            except:
                warnings.warn(
                    'Error when estimating decaying zerocrossing after peak ' +
                    str(i) + ' at sample ' + str(Ps[i]) +
                    '. Therefore, the zerocrossing has been set to halfway between the two extrema.'
                )
                zeroxD[i] = Ps[i] + int(len(x_temp) / 2.)

    return zeroxR, zeroxD