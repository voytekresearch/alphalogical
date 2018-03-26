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
cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

# --------------------------------------------------------


def find_closest_idx(t, times):
    times = np.asarray(times)
    idx = (np.abs(times - t)).argmin()

    return idx


def find_closest(t, times):
    idx = find_closest_idx(t, times)
    return times[idx]


@memory.cache
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


def find_smith_bursts(x_ref, x_lfp, center, fs):
    # Power
    alpha = butter_bandpass_filter(x_lfp, center - 2, center + 2, fs, order=2)
    alpha_pow = np.abs(hilbert(alpha))

    # Est. thresh
    # Find low
    low_range = [center - 4, center - 2]
    if low_range[0] < 0:
        low_range[0] = 1
    if low_range[1] < 2:
        low_range[1] = 2
    low = np.abs(
        hilbert(
            butter_bandpass_filter(
                x_ref, low_range[0], low_range[1], fs, order=2)))

    # Find high
    high_range = [center + 2, center + 4]
    high = np.abs(
        hilbert(
            butter_bandpass_filter(
                x_ref, high_range[0], high_range[1], fs, order=2)))

    # Set thresh
    M = np.mean((low + high) / 2.0)

    # Find bursts
    bursts = find_bursts(alpha_pow, 3 * M, 1.5 * M)

    return bursts, alpha, alpha_pow, M


def analyze_smith_bursts(results_code, bursts, alpha_pow, hg_pow, rates, fs):

    # Init
    n_bursts = len(bursts)
    results = defaultdict(list)

    for k, b in enumerate(bursts):
        burst_time = float(len(b)) / fs
        alpha_power = np.median(alpha_pow[b])
        hg_power = np.median(hg_pow[b])
        rate = np.median(rates[b])

        # Find peak/troughs
        #     i, j = b.min(), b.max()

        #     peaks, troughs = find_extrema(
        #         alpha[i:j], fs, 
        #         (center - bw, center + bw), 
        #         boundary=None, first_extrema='peak'
        #     )

        #     n_cycles = len(peaks)

        # Peak mask
        # TODO stats

        # Trough mask
        # TODO stats

        # Save results       
        results["results_code"].append(results_code)
        results["burst_index"].append(k)
        results["n_bursts"].append(n_bursts)
        results["burst_time"].append(burst_time)
        results["alpha_power"].append(alpha_power)
        results["hg_power"].append(hg_power)
        results["rate"].append(rate)
    #     results["n_cycles"].append(n_cycles)
    #     results["peak_alpha_power"].append(peak_alpha_power)
    #     results["peak_hg_power"].append(peak_hg_power)
    #     results["peak_rate"].append(peak_rate)
    #     results["trough_alpha_power"].append(trough_alpha_power)
    #     results["trough_hg_power"].append(trough_hg_power)
    #     results["trough_rate"].append(trough_rate)                                

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
        Ps[p] = np.argmax(x[mrzerorise:nfzerofall]) + mrzerorise

    # Calculate trough samples
    Ts = np.zeros(T, dtype=int)
    for tr in range(T):
        # Calculate the sample range between the most recent zero fall
        # and the next zero rise
        mrzerofall = zerofallN[tr]
        nfzerorise = zeroriseN[zeroriseN > mrzerofall][0]
        # Identify time of trough
        Ts[tr] = np.argmin(x[mrzerofall:nfzerorise]) + mrzerofall

    # Remove peaks and troughs within the boundary limit
    Ps = Ps[np.logical_and(Ps > boundary, Ps < len(x) - boundary)]
    Ts = Ts[np.logical_and(Ts > boundary, Ts < len(x) - boundary)]

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