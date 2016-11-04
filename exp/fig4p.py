# -*- coding: utf8 -*-
"""Usage: fig4p.py NAME STIM
    [-t T] 
    [-n N] 
    [-l L] 

Phase effects with stimulus presentation

    Arguments:
        NAME        name of the result file
        STIM        config file with stimulus presentation

    Options:
        -h --help    show this screen
        -t T         fraction of std dev for detection [default: 1]
        -n N         number of trials [default: 360]
        -l L         the number of power levels [default: 6]

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
if __name__ == "__main__":
    args = docopt(__doc__, version='alpha')
                                                  
    name = args['NAME']

    base_par = args['STIM']
    par_path = os.path.dirname(base_par)
    par_name = os.path.splitext(os.path.basename(base_par))[0]

    threshold = float(args['-t'])

    # -------------------------------------------------------------------
    t = 0.8
    dt = 1e-3
    times = create_times(t, dt)

    n_stim = int(args['-n'])
    # rate = 20 

    pow1 = 1

    n_pow2 = int(args['-l'])
    powers2 = np.linspace(1, 6, n_pow2)

    sigma = 10  
    loc = ['r_E', ]

    # Init two sets of detect data variables.
    # 1. suppression model.
    hits = np.zeros((n_stim, len(powers2)))
    false_alarms = np.zeros((n_stim, len(powers2)))
    d_primes = np.zeros((n_stim, len(powers2)))

    # 2. gatin' model
    p_lefts = np.zeros((n_stim, len(powers2)))
    p_rights = np.zeros((n_stim, len(powers2)))

    for i, pow2 in enumerate(powers2):
        # Load params for first power level
        n_pop, pops, names, inputs, backs, conns = parse(base_par)

        # Override background rate
        backs[0]['re'] *= pow1
        par1 = os.path.join(par_path, "{}_pow1{}.yaml".format(par_name, pow2))
        save(par1, n_pop, pops, inputs, backs, conns)

        # Reload params for power level 2
        n_pop, pops, names, inputs, backs, conns = parse(base_par)

        # Override background rate
        backs[0]['re'] *= pow2
        par2 = os.path.join(par_path, "{}_pow2{}.yaml".format(par_name, pow2))
        save(par2, n_pop, pops, inputs, backs, conns)
        
        # No stim
        inputs[0]['mode']['args'][0] = 0.00001
        par3 = os.path.join(par_path, "{}_pow2{}_nostim.yaml".format(
            par_name, pow2))
        save(par3, n_pop, pops, inputs, backs, conns)  

        # Trial loop
        for j in range(n_stim):
            print(">>> Name {} (stim), power level 2 {}, stim {}".format(
                name, pow2, j))

            # ---
            # integrate
            print(">>> Power level 1")
            ys1, layers1, phi1, rate1, stim1, params1 = run(None,
                                                      times,
                                                      par1,
                                                      sigma=sigma,
                                                      loc=loc,
                                                      stim_seed=j)
           # integrate
            print(">>> Power level 2")
            ys2, layers2, phi2, rate2, stim2, params2 = run(None,
                                                      times,
                                                      par2,
                                                      sigma=sigma,
                                                      loc=loc,
                                                      stim_seed=j)
            # integrate
            print(">>> Power level 2, no stim")
            ys3, layers3, phi3, rate3, stim3, params3 = run(None,
                                                      times,
                                                      par3,
                                                      sigma=sigma,
                                                      loc=loc,
                                                      stim_seed=j)

            idx1 = params1.idx
            idx2 = params2.idx
            idx3 = params3.idx

            # E activity
            re1 = ys1[:, idx1['r_E']]
            re2 = ys2[:, idx2['r_E']]
            re3 = ys3[:, idx3['r_E']]

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

            M_pow1 = re1[m_post].mean()
            M_pow2 = re2[m_post].mean()
            M_pow3 = re3[m_post].mean()

            SD_pow1 = re1[m_post].std()
            SD_pow2 = re2[m_post].std()
            SD_pow3 = re3[m_post].std()

            # -
            # Abs detection with decreasing gain 
            # compare pow2 (with stim) to
            # pow3. which has no stim
            #
            # Must be one of the two below.
            hit = 0
            false_alarm = 0

            # Hit?
            if M_pow2 > (M_pre + (SD_pre * threshold)):
                hit = 1

            # False alarm?
            if M_pow3 > (M_pre + (SD_pre * threshold)):
                false_alarm = 1

            hits[j, i] = hit
            false_alarms[j, i] = false_alarm
            
            # D prime
            d_prime = (M_pow2 - M_pow3) / np.sqrt((SD_pow2**2 + SD_pow3**2) / 2)
            d_primes[j, i] = d_prime

            # --
            # diff detection (pow1 versus pow2)
            # calculate the p(left) ane p(right)

            # pow1 := left
            # pow2 := right
            M_left = M_pow1
            M_right = M_pow2

            p_left = np.exp(M_left) / np.sum([np.exp(M_left), np.exp(M_right)])
            p_right = 1 - p_left
            
            p_lefts[j, i] = p_left
            p_rights[j, i] = p_right

    # -------------------------------------------------------------------
    # Calculate remaining signal detection statistics
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
        p_lefts=p_lefts,
        p_rights=p_rights,
        d_primes=d_primes,
        n_stim=n_stim,
        stims=range(n_stim),
        pow1=pow1,
        powers2=powers2,
        sigma=sigma,
        loc=loc)
