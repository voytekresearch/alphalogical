from __future__ import division
import os
import sys
import numpy as np
import pyentropy as en
from bluemass.bm import run
from fakespikes.util import create_times
from pykdf.kdf import save_kdf


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
def norm(x):
    return (x - x.min()) / (x.max() - x.min())


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
save_path = sys.argv[1]
pars_path = sys.argv[2]
exp = sys.argv[3]

t = 2
dt = 1e-3
drop_before = 0.5
times = create_times(t, dt)

n_stim = 36
loc = [
    'r_E',
]
m = 10  # 8 is Ince's advice

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# a + b. information flow and PAC
n_sigma = 10
sigmas = np.linspace(1, 50, n_sigma)

# -------------------------------------------------------------------
# a. async
if exp == "a":
    name = "a"

    nopac = os.path.join(pars_path, 'pac_async_noalpha.yaml')
    pac = os.path.join(pars_path, 'pac_async.yaml')

    del_mi = np.zeros((n_stim, n_sigma))
    for i, s in enumerate(sigmas):
        for j in range(n_stim):
            print(">>> Fig {}, sigma {}, stim {}".format(name, s, j))

            # ---
            # integrate
            ys, layers, phi, rate, stim, params = run(None,
                                                      times,
                                                      nopac,
                                                      sigma=s,
                                                      loc=loc,
                                                      stim_seed=j)

            idx = params.idx

            # extract result
            x_in = np.asarray([stim(t, 0) for t in times])
            x_e = ys[:, idx['r_E']].flatten()

            # drop burn in time
            x_in = x_in[times > drop_before]
            x_e = x_e[times > drop_before]

            # normalize
            x_in = norm(x_in)
            x_e = norm(x_e)

            # prep MI variables
            to_calc = ('HX', 'HY', 'HXY')
            q_in, _, _ = en.quantise(x_in, m)
            q_e, _, _ = en.quantise(x_e, m)

            # MI
            info = en.DiscreteSystem(q_in, (1, m), q_e, (1, m))
            info.calculate_entropies(method='pt', calc=to_calc)
            mi_no = info.I()

            # ---
            ys, layers, phi, rate, stim, params = run(None,
                                                      times,
                                                      pac,
                                                      sigma=s,
                                                      loc=loc,
                                                      stim_seed=j)

            idx = params.idx

            x_in = np.asarray([stim(t, 0) for t in times])
            x_e = ys[:, idx['r_E']].flatten()

            x_in = x_in[times > drop_before]
            x_e = x_e[times > drop_before]

            to_calc = ('HX', 'HY', 'HXY')
            q_in, _, _ = en.quantise(x_in, m)
            q_e, _, _ = en.quantise(x_e, m)

            info = en.DiscreteSystem(q_in, (1, m), q_e, (1, m))
            info.calculate_entropies(method='pt', calc=to_calc)
            mi_pac = info.I()

            # Save
            del_mi[j, i] = mi_pac - mi_no

    save_kdf(
        os.path.join(save_path, name),
        del_mi=del_mi,
        n_stim=n_stim,
        n_sigma=n_sigma,
        stims=range(n_stim),
        sigmas=sigmas)

# -------------------------------------------------------------------
# b. osc
elif exp == "b":
    name = "b"
    nopac = os.path.join(pars_path, 'pac_osc_noalpha.yaml')
    pac = os.path.join(pars_path, 'pac_osc.yaml')

    del_mi = np.zeros((n_stim, n_sigma))
    for i, s in enumerate(sigmas):
        for j in range(n_stim):
            print(">>> Fig {}, sigma {}, stim {}".format(name, s, j))

            # ---
            # integrate
            ys, layers, phi, rate, stim, params = run(None,
                                                      times,
                                                      nopac,
                                                      sigma=s,
                                                      loc=loc,
                                                      stim_seed=j)

            idx = params.idx

            # extract result
            x_in = np.asarray([stim(t, 0) for t in times])
            x_e = ys[:, idx['r_E']].flatten()

            # drop burn in time
            x_in = x_in[times > drop_before]
            x_e = x_e[times > drop_before]

            # normalize
            x_in = norm(x_in)
            x_e = norm(x_e)

            # prep MI variables
            to_calc = ('HX', 'HY', 'HXY')
            q_in, _, _ = en.quantise(x_in, m)
            q_e, _, _ = en.quantise(x_e, m)

            # MI
            info = en.DiscreteSystem(q_in, (1, m), q_e, (1, m))
            info.calculate_entropies(method='pt', calc=to_calc)
            mi_no = info.I()

            # ---
            ys, layers, phi, rate, stim, params = run(None,
                                                      times,
                                                      pac,
                                                      sigma=s,
                                                      loc=loc,
                                                      stim_seed=j)

            idx = params.idx

            x_in = np.asarray([stim(t, 0) for t in times])
            x_e = ys[:, idx['r_E']].flatten()

            x_in = x_in[times > drop_before]
            x_e = x_e[times > drop_before]

            to_calc = ('HX', 'HY', 'HXY')
            q_in, _, _ = en.quantise(x_in, m)
            q_e, _, _ = en.quantise(x_e, m)

            info = en.DiscreteSystem(q_in, (1, m), q_e, (1, m))
            info.calculate_entropies(method='pt', calc=to_calc)
            mi_pac = info.I()

            # Save
            del_mi[j, i] = mi_pac - mi_no

    save_kdf(
        os.path.join(save_path, name),
        del_mi=del_mi,
        n_stim=n_stim,
        n_sigma=n_sigma,
        stims=range(n_stim),
        sigmas=sigmas)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# c. phase diff
# two populations, calculate the MI between their async E as a function
# of phase
# if gain frequency is 10 Hz (and it is), a temporal offset of 50 ms
# will take two 10 Hz oscillators 180 degrees out of phase.
# while a 100 ms lag will realign them
elif exp == "c":
    name = "c"

    sigma = 10

    phases = [0, 5, 11, 22, 45, 90, 180]  # degrees, obvs
    n_phase = len(phases)
    par_phases = [
        'phase_diff_ph0.yaml', 'phase_diff_ph5.yaml', 'phase_diff_ph11.yaml',
        'phase_diff_ph22.yaml', 'phase_diff_ph45.yaml', 'phase_diff_ph90.yaml',
        'phase_diff_ph180.yaml'
    ]

    mi_phase = np.zeros((n_stim, n_phase))
    for i in range(n_stim):
        # Create a reference standard for rE
        ys, layers, phi, rate, stim, params = run(
            None,
            times,
            os.path.join(pars_path, 'phase_diff_ph0.yaml'),
            sigma=sigma,
            loc=loc,
            stim_seed=i)

        idx = params.idx

        x_std = ys[:, idx['r_E']].flatten()
        x_std = x_std[times > drop_before]

        # normalize
        x_std = norm(x_std)

        for j, (phase, par) in enumerate(zip(phases, par_phases)):
            print(">>> Fig {}, stim {}, phase {}".format(name, i, phase))

            ys, layers, phi, rate, stim, params = run(
                None,
                times,
                os.path.join(pars_path, par),
                sigma=sigma,
                loc=loc,
                stim_seed=i)

            idx = params.idx

            x_e = ys[:, idx['r_E']].flatten()
            x_e = x_e[times > drop_before]

            # normalize
            x_e = norm(x_e)

            # prep MI variables
            to_calc = ('HX', 'HY', 'HXY')
            q_std, _, _ = en.quantise(x_std, m)
            q_e, _, _ = en.quantise(x_e, m)

            # MI
            info = en.DiscreteSystem(q_std, (1, m), q_e, (1, m))
            info.calculate_entropies(method='pt', calc=to_calc)

            mi_phase[i, j] = info.I()

    save_kdf(
        os.path.join(save_path, name),
        sigma=sigma,
        mi_phase=mi_phase,
        n_stim=n_stim,
        n_phase=n_phase,
        phases=phases,
        par_phases=par_phases)

else:
    raise ValueError("exp must be a, b or c")
