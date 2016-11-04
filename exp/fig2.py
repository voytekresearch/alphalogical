import sys, os
import numpy as np
from bluemass.fi import lif
from pykdf.kdf import save_kdf, load_kdf
import fakespikes.util as sp

save_path = sys.argv[1]

# Const params
t = 1
n_trials = 30

Is = np.linspace(0, 50e-3, 100)

r_e = 135
r_i = 135
min_rate = 30

w_e = 4e-9
w_i = 3.91 * w_e
f = 10

# --
# a)
print("Running a)")
amps = [1, 2, 3]

n_bursts = None

fis = []
for a in amps:
    print("Amp. factor {}".format(a))

    fi_t = []
    for k in range(n_trials):
        fi_t.append(
            lif(t,
                Is,
                f,
                r_e=r_e * a,
                r_i=r_i * a,
                w_e=w_e,
                w_i=w_i,
                min_rate=min_rate * a, # TODO keep this scaling too?
                n_bursts=n_bursts,
                back_seed=42 * k,
                verbose=False))

    fis.append(np.vstack(fi_t).mean(0))

save_kdf(
    os.path.join(save_path, "a"), Is=Is, f=f, amps=amps, fis=fis)

# --
# b)
print("Running b)")
n_bursts = [2, 4, 6, 8]
fis = []
for n in n_bursts:
    print("N bursts {}".format(n))

    fi_t = []
    for k in range(n_trials):
        fi_t.append(
            lif(t,
                Is,
                f,
                r_e=r_e,
                r_i=r_i,
                w_e=w_e,
                w_i=w_i,
                min_rate=min_rate,
                n_bursts=n,
                back_seed=42 * k,
                verbose=False))

    fis.append(np.vstack(fi_t).mean(0))

save_kdf(os.path.join(save_path, "b"), Is=Is, f=f, n_bursts=n_bursts, fis=fis)
