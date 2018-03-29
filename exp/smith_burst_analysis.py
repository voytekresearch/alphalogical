import fire
import os
import h5py
import numpy as np

from alphalogical.util import butter_bandpass_filter
from alphalogical.util import hilbert_power

from alphalogical.util import load_good_smith_files
from alphalogical.util import load_smith_foof_results
from alphalogical.util import load_smith_rates
from alphalogical.util import load_smith_lfps

from alphalogical.util import find_bursts
from alphalogical.util import find_extrema

from alphalogical.util import analyze_bursts
from alphalogical.util import shift_bursts
from alphalogical.util import shift_extrema
from alphalogical.util import analyze_extrema
from alphalogical.util import write_analysis


def run(name, data_path, n=None, percent_segment=1, bw=2, verbose=False):
    # --------------------------------------------------------------
    if verbose:
        print(">>> Finding Smith files".format(name))
    lfp_files, rate_files, seg_files = load_good_smith_files(data_path)

    # Run only file n?
    if n is not None:
        lfp_files = [lfp_files[n]]
        rate_files = [rate_files[n]]
        sef_files = [seg_files[n]]

    # --------------------------------------------------------------
    # For each data file, get the segments.
    for n_file, (lfp_file, rate_file,
                 seg_file) in enumerate(zip(lfp_files, rate_files, seg_files)):

        # Load
        header, segments = load_smith_foof_results(seg_file)
        if verbose:
            print(">>> Loading lfp file {}".format(lfp_file))
            print(">>> Loading rate file {}".format(rate_file))
            print(">>> Loading segment file {}".format(seg_file))
            print(">>> Segment header: {}".format(header))
            print(">>> Segment number: {}".format(len(segments)))
            print(">>> Segment example row: {}".format(segments[0]))
            print(">>> Running segment analysis....")

        if percent_segment < 1:
            n_segments = int(len(segments) * percent_segment)
            if verbose:
                print(">>> Sampling {}/{} segments".format(n_segments,
                                                           len(segments)))
            segments = segments[0:n_segments]

        # --------------------------------------------------------------
        # Analyze each segment for bursts/peak/troughs
        for seg in segments:
            m, c, i, j, center, power, _ = seg

            if verbose:
                print(">>> Loading data for segment {}.".format(m))

            fs, times, lfp = load_smith_lfps(lfp_file, c)
            T = times.max()

            # Load rate
            _, rate = load_smith_rates(rate_file, c, (0, T), 1 / fs)

            # pad rate to match lfp (1 bin correction)
            rate = np.pad(rate, (1, 0), 'constant')

            assert lfp.shape == rate.shape, "lfp and rate don't match"

            # --------------------------------------------------------------
            if verbose:
                print(">>> Analyzing segment {}.".format(m))

            # Do filtering over the whole electrode, for consistency
            # -
            # Alpha
            alpha_range = (center - bw, center + bw)
            alpha = butter_bandpass_filter(
                lfp, alpha_range[0], alpha_range[1], fs, order=2)

            # Watch out for even signals....
            if alpha.size % 2 == 0:
                alpha_pow = hilbert_power(alpha)
            else:
                # Make even
                alpha_pow = hilbert_power(alpha[1:])
                # then pad
                alpha_pow = np.pad(alpha_pow, (0, 1), 'constant')

            # -
            # HG
            hg = butter_bandpass_filter(lfp, 100, 150, fs, order=2)

            # Watch out for even signals....
            if hg.size % 2 == 0:
                hg_pow = hilbert_power(hg)
            else:
                # Make even
                hg_pow = hilbert_power(hg[1:])
                # then pad
                hg_pow = np.pad(hg_pow, (0, 1), 'constant')

            hg_pow = hilbert_power(hg)

            assert alpha.shape == lfp.shape, "lfp and alpha don't match"
            assert hg.shape == lfp.shape, "lfp and HG don't match"

            # --------------------------------------------------------------
            # Find low
            low_range = [center - (bw + 2), center - bw]
            if low_range[0] < 0:
                low_range[0] = 1
            if low_range[1] < 2:
                low_range[1] = 2
            low = butter_bandpass_filter(
                lfp, low_range[0], low_range[1], fs, order=2)
            low = hilbert_power(low)

            # Find high
            high_range = [center + bw, center + (bw + 2)]
            high = butter_bandpass_filter(
                lfp, high_range[0], high_range[1], fs, order=2)
            high = hilbert_power(high)

            # --------------------------------------------------------------
            # Only run for this mth segment...
            # Bursts
            M = np.mean((low + high) / 2.0)
            bursts = find_bursts(alpha_pow[i:j], 3 * M, 1.5 * M)

            # Peaks
            peaks, troughs = find_extrema(
                alpha[i:j],
                fs,
                alpha_range,
                boundary=None,
                first_extrema='peak')

            # Shift indices to match data
            bursts = shift_bursts(bursts, i)
            peaks = shift_extrema(peaks, i)
            troughs = shift_extrema(troughs, i)

            # --------------------------------------------------------------
            if verbose:
                print(">>> Saving results for segment {}.".format(m))

            # alpha
            results = analyze_bursts("{}_{}".format(n_file, m), bursts,
                                     alpha_pow)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_alpha_burst_results"),
                    results,
                    mode="a")

            results = analyze_extrema(
                "{}_{}".format(n_file, m), peaks, alpha_pow, n_sample=10)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_alpha_peak_results"),
                    results,
                    mode="a")

            results = analyze_extrema(
                "{}_{}".format(n_file, m), troughs, alpha_pow, n_sample=10)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_alpha_trough_results"),
                    results,
                    mode="a")

            # hg
            results = analyze_bursts("{}_{}".format(n_file, m), bursts, hg_pow)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_hg_burst_results"),
                    results,
                    mode="a")

            results = analyze_extrema(
                "{}_{}".format(n_file, m), peaks, hg_pow, n_sample=10)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_hg_peak_results"),
                    results,
                    mode="a")

            results = analyze_extrema(
                "{}_{}".format(n_file, m), troughs, hg_pow, n_sample=10)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_hg_trough_results"),
                    results,
                    mode="a")

            # rate
            results = analyze_bursts("{}_{}".format(n_file, m), bursts, rate)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_rate_burst_results"),
                    results,
                    mode="a")

            results = analyze_extrema(
                "{}_{}".format(n_file, m), peaks, rate, n_sample=10)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_rate_peak_results"),
                    results,
                    mode="a")

            results = analyze_extrema(
                "{}_{}".format(n_file, m), troughs, rate, n_sample=10)
            if results is not None:
                write_analysis(
                    os.path.join(data_path, name + "_rate_trough_results"),
                    results,
                    mode="a")

    return None


if __name__ == "__main__":
    fire.Fire(run)
