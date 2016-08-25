from tanh_neuron import TanhWithBias

import numpy as np
import matplotlib.pyplot as plt
import ipdb

SEED = 0
nrn = TanhWithBias(seed=SEED)
n_neurons = 5
step_size = 0.1
step_list = list(np.arange(-1, 1, step_size))
steps = len(step_list)
res = np.zeros((steps, n_neurons))

for step, x in enumerate(step_list):
    new_val = np.zeros(n_neurons)
    in_j = x * np.ones(n_neurons)
    nrn.step_math(0.001, in_j, new_val)
    res[step, :] = new_val

plt.plot(res)
plt.show()
ipdb.set_trace()
