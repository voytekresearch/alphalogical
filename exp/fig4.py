# -*- coding: utf8 -*-
"""Usage: fig4.py NAME STIM
    [-t THRESHOLD] 

Phase effects with stimulus presentation

    Arguments:
        NAME        name of the result file
        STIM        config file with stimulus presentation

    Options:
        -h --help               show this screen
        -t T                    fraction of std dev for detection [default: 1]

"""
from __future__ import division
import os
import sys
import numpy as np
import pyentropy as en
from docopt import docopt
from bluemass.bm import run
from bluemass.params import parse, save
from fakespikes.util import create_times
from pykdf.kdf import save_kdf
from scipy.stats import norm


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
if __name__ == "__main__":
    args = docopt(__doc__, version='alpha')
                                                  
    name = args['NAME']

    base_par1 = args['STIM']
    par1_path = os.path.dirname(base_par1)
    par1_name = os.path.splitext(os.path.basename(base_par1))[0]

    threshold = float(args['-t'])

    # -------------------------------------------------------------------
    t = 0.8
    dt = 1e-3
    times = create_times(t, dt)

    n_stim = 36 * 10 # 36 trials 10 per block in Mathewson 
    rates = np.linspace(1, 100, 10) 

    sigma = 10  
    loc = ['r_E', ]

    hits = np.zeros((n_stim, len(rates)))
    false_alarms = np.zeros((n_stim, len(rates)))
    d_primes = np.zeros((n_stim, len(rates)))
    for i, r in enumerate(rates):
        
        # Load params
        n_pop, pops, names, inputs, backs, conns = parse(base_par1)

        # Override stim
        inputs[0]['mode']['args'][0] = r
        par1 = os.path.join(par1_path, "{}_r{}.yaml".format(par1_name, r))
        save(par1, n_pop, pops, inputs, backs, conns)

        # No stim
        inputs[0]['mode']['args'][0] = 0.00001
        par2 = os.path.join(par1_path, "{}_r{}_nostim.yaml".format(par1_name, r))
        save(par2, n_pop, pops, inputs, backs, conns)    

        for j in range(n_stim):

            print(">>> Name {} (stim), rate {}, stim {}".format(name, r, j))

            # ---
            # integrate
            print(">>> Stimulus")
            ys1, layers1, phi1, rate1, stim1, params1 = run(None,
                                                      times,
                                                      par1,
                                                      sigma=sigma,
                                                      loc=loc,
                                                      stim_seed=j)
           # integrate
            print(">>> No stimulus")
            ys2, layers2, phi2, rate2, stim2, params2 = run(None,
                                                      times,
                                                      par2,
                                                      sigma=sigma,
                                                      loc=loc,
                                                      stim_seed=j)

            idx1 = params1.idx
            idx2 = params2.idx

            # E activity
            re1 = ys1[:, idx1['r_E']]
            re2 = ys2[:, idx1['r_E']]

            # After skiping the fiest 100 ms as burn in
            # locate activity for the 200 ms preceeding
            # stimulus onset. This is the prestim
            # comparison period.
            m_pre = np.logical_and(times > 0.1, times < 0.3)  

            # stim onset at 0.4 s; we take 20 ms before and after the stim
            # as the detection windows.
            m_post = np.logical_and(times > 0.4, times < 0.42)  

            M_pre = np.mean([re1[m_pre].mean(), re2[m_pre].mean()])
            SD_pre = np.mean([re1[m_pre].std(), re2[m_pre].std()])

            M_stim = re1[m_post].mean()
            M_no = re2[m_post].mean()

            SD_stim = re1[m_post].std()
            SD_no = re2[m_post].std()

            # Detect this trial?
            # Must be one of the two below.
            hit = 0
            false_alarm = 0

            # Hit?
            if M_stim > (M_pre + (SD_pre * threshold)):
                hit = 1

            # False alarm?
            if M_no > (M_pre + (SD_pre * threshold)):
                false_alarm = 1

            hits[j, i] = hit
            false_alarms[j, i] = false_alarm
            
            # D prime
            d_prime = (M_stim - M_no) / np.sqrt((SD_stim**2 + SD_no**2) / 2)
            d_primes[j, i] = d_prime

    # -------------------------------------------------------------------
    # Calculate signal detection statistics

    # First calculate misses and correct_rejections
    misses = np.ones_like(hits) - hits
    correct_rejections = np.ones_like(false_alarm) - false_alarms
 
    # -------------------------------------------------------------------
    # Save
    save_kdf(
        name,
        hits=hits,
        misses=misses,
        correct_rejections=correct_rejections,
        false_alarms=false_alarms,
        d_primes=d_primes,
        n_stim=n_stim,
        stims=range(n_stim),
        rates=rates)

