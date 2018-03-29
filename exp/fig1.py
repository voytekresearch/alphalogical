import sys, os
import numpy as np
from bluemass.fi import lif
from pykdf.kdf import save_kdf, load_kdf
import fakespikes.util as sp

# --
save_path = sys.argv[1]

# Const params
t = 60

Is = np.linspace(0, 50e-3, 100)
w_e = 4e-9
w_i = 3.91 * w_e

# --
# a) on/off
r_e = 135
r_i = 135
min_rate = 30

f = 0
fi1, trains1, ge1, gi1, v1 = lif(t,
                                 Is,
                                 f,
                                 r_e=r_e,
                                 r_i=r_i,
                                 w_e=w_e,
                                 w_i=w_i,
                                 min_rate=min_rate,
                                 return_trains=True)
ns1, ts1 = sp.spikedict_to(trains1)

f = 10
fi2, trains2, ge2, gi2, v2 = lif(t,
                                 Is,
                                 f,
                                 r_e=r_e,
                                 r_i=r_i,
                                 w_e=w_e,
                                 w_i=w_i,
                                 min_rate=min_rate,
                                 return_trains=True)
ns2, ts2 = sp.spikedict_to(trains2)

save_kdf(
    os.path.join(save_path, "a"),
    Is=Is,
    f=10,
    fixed=1,
    osc=2,
    fi1=fi1,
    fi2=fi2,
    ns1=ns1,
    ts1=ts1,
    ns2=ns2,
    ts2=ts2,
    ge1=ge1,
    gi1=gi1,
    v1=v1,
    ge2=ge2,
    gi2=gi2,
    v2=v2, )

# --
# b) phase plot, and full FI for peak/trough
# we turn osc off because were estimating gain at peak/trough by
# simulating at the peak/trough rates for longer periods of time than 
# the peak/trough lasts peak.
f = 0

rs = np.linspace(30, 135, 6)
phases = np.linspace(0, 180, 6)

fi_phase = []
for r in rs:
    r_e = r
    r_i = r
    fi_r, _, _, _, _ = lif(t,
                           Is,
                           f,
                           r_e=r_e,
                           r_i=r_i,
                           w_e=w_e,
                           w_i=w_i,
                           min_rate=min_rate,
                           return_trains=True)

    fi_phase.append(fi_r.mean())
fi_phase = np.asarray(fi_phase)

# peak
r_e = 135
r_i = 135
fi1, trains1, ge1, gi1, v1 = lif(t,
                                 Is,
                                 f,
                                 r_e=r_e,
                                 r_i=r_i,
                                 w_e=w_e,
                                 w_i=w_i,
                                 min_rate=min_rate,
                                 return_trains=True)
ns1, ts1 = sp.spikedict_to(trains1)

# trough
r_e = 30
r_i = 30
fi2, trains2, ge2, gi2, v2 = lif(t,
                                 Is,
                                 f,
                                 r_e=r_e,
                                 r_i=r_i,
                                 w_e=w_e,
                                 w_i=w_i,
                                 min_rate=min_rate,
                                 return_trains=True)
ns2, ts2 = sp.spikedict_to(trains2)

save_kdf(
    os.path.join(save_path, "b"),
    phases=phases,
    fi_phase=fi_phase,
    peak=1,
    trough=2,
    Is=Is,
    fi1=fi1,
    fi2=fi2,
    ns1=ns1,
    ts1=ts1,
    ns2=ns2,
    ts2=ts2,
    ge1=ge1,
    gi1=gi1,
    v1=v1,
    ge2=ge2,
    gi2=gi2,
    v2=v2)
