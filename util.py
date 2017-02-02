from copy import deepcopy
import numpy as np
from scipy.signal import butter, lfilter


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
            burst = [i, ]
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


if __name__ == "__main__":
    from scipy.io import loadmat
    from scipy.signal import hilbert
    from scipy.signal import butter, lfilter
    from alphalogical.util import *

    fs = 1000
    t = 120  # seconds long
    alpha = loadmat('data/alpha_data.mat')
    times = np.linspace(0, 120, 120 * 500)

    task = alpha['oz_dat_task'][0, :]

    # smooth a little

    data = butter_bandpass_filter(task, .1, 30, fs, order=2)

    data_f = butter_bandpass_filter(data, 8, 12, fs, order=2)
    data_pow = np.abs(hilbert(data_f))

    low = np.abs(hilbert(butter_bandpass_filter(data, 5, 7, fs, order=2)))
    high = np.abs(hilbert(butter_bandpass_filter(data, 14, 16, fs, order=2)))
    base = (low + high) / 2.0

    task = butter_bandpass_filter(task, .1, 30, 500, order=2)

    # Find burst indices
    start = 1.5
    stop = 1.0
    task_b = find_bursts(data_pow, base.mean() * start, base.mean() * stop)

    # Burst stats
    print(len(task_b))