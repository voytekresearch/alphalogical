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


def find_bursts2(x, burst_thresh, zero_thresh):
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
    for i in range(1, n):
        z_last = zeros[i - 1]
        c = candidates[i]

        # Still bursting?
        if b and c:
            burst.append(i)
        # Burst onset?
        elif c and z_last:
            burst = [i, ]
            b = True
        # Burst over.
        else:
            # Store that bursts index
            bursts.append(deepcopy(burst))

            # and reset
            burst = []
            b = False

    return bursts
