from utils import get_conceptors

import scipy.io
import numpy as np
import ipdb


sig_num = 2
t_steps = 451
apert = np.ones(sig_num) * 10
n_neurons = 600

patts = []
concs = []

for fi in range(1, sig_num):
    patts.append(scipy.io.loadmat("test_files/patt%s.mat" % fi)['xCollector'].T)
    concs.append(scipy.io.loadmat("test_files/conc%s.mat" % fi)['af'])

corrs = scipy.io.loadmat("test_files/corrs.mat")['patternRs'][0]
cv = scipy.io.loadmat("test_files/conc_vals.mat")['Cs']

res_c = get_conceptors(patts, sig_num-1, t_steps, apert, n_neurons)

ipdb.set_trace()
