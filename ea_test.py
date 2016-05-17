# test architecture with reduced functions

import scipy.io
from scipy import interpolate
import nengo
import numpy as np
from dmp_utils import *
import matplotlib.pyplot as plt
import ipdb

def make_dmp_net(functions, input_obj, output_obj, name=""):
    """create one point attractor per dimension with goals as nodes
    and one unified neuron input for inhibition

    TODO:
        - make the connections for inhibition (proxy for output of BG)
        - make the connections for scaling (proxy for output of Thal)
    """
    n = nengo.Network(label=name)
    n.pt_attractors = []
    n.conn_funcs = []
    n.f_conns = []
    with n:
        n.output = nengo.Node(size_in=len(functions))

        for d in range(len(functions)):
            goal = nengo.Node([0], label="goal_%s" %(d))
            attractor = gen_point_attractor(n, goal, n_neurons=500)
            attractor.label = "pt_attr_%s" %(d)
            nengo.Connection(attractor[0], n.output[d], synapse=None)
            n.pt_attractors.append(attractor)

    for f_i, func in enumerate(functions):
        dest = func(np.linspace(-np.pi, np.pi, ea_func_steps))
        dest = dest.reshape((-1, 1))
        force_func = gen_forcing_functions(dest)[0]
        n.conn_funcs.append(lambda x, force_func=force_func: force(x, force_func))

        # there's still the little bump, but it doesn't seem as bad?
        # it's just that the force can't change instanteously
        # it seems like there should be some way to get around that
        # or maybe it won't matter for large force?
        n.f_conns.append(nengo.Connection(input_obj, n.pt_attractors[f_i][1], synapse=None,
                         function=n.conn_funcs[f_i]))

    nengo.Connection(n.output, output_obj, synapse=None)
    return n


pattern_file_names = [
     "nnRawExaStride",
     "nnRawSlowWalk",
     "nnRawWalk",
     "nnRawRunJog",
     "nnRawCartWheel",
     "nnRawWaltz",
     "nnRawCrawl",
     "nnRawStandup",
     "nnRawGetdown",
     "nnRawSitting",
     "nnRawGetSeated",
     "nnRawStandupFromStool",
     "nnRawBox1",
     "nnRawBox2",
     "nnRawBox3",
]

output_dims = 2
pattern_num = 1
pattern_file_names = pattern_file_names[:pattern_num]

function_list = [
    [
        lambda x: np.sin(x),
        lambda x: -np.sin(x)
    ]
]

ea_n_neurons = 300
ea_func_steps = 100
np.random.seed(3)
# maps from input value (in this case, theta) to output value

model = nengo.Network()
tau = 0.1
dmp_net_list = []

with model:
    osc = nengo.Ensemble(n_neurons=1, dimensions=3, neuron_type=nengo.Direct())

    def cycle(x):
        """makes a speed controlled oscillator"""
        a = 1.0
        b = 2.0 * np.pi * x[2]
        r = np.sqrt(x[0]**2.0 + x[1]**2.0)
        theta = np.arctan2(x[1], x[0])
        dr = 10.0*(-r**3.0 + a*r)
        dtheta = b

        dx = dr*np.cos(theta) - r*np.sin(theta)*dtheta
        dy = dr*np.sin(theta) + r*np.cos(theta)*dtheta

        return [x[0] + tau*dx, x[1] + tau*dy]

    nengo.Connection(osc, osc[:2], synapse=tau, function=cycle)

    rate = nengo.Node([1])
    nengo.Connection(rate, osc[2])

    bump = nengo.Node(lambda t: 1 if t < 0.5 else 0)
    nengo.Connection(bump, osc[0])

    # controllers
    #inhibit_control = nengo.Node([0])
    #scale_control = nengo.Node([0]*pattern_num)

    output = nengo.networks.EnsembleArray(n_neurons=1, n_ensembles=output_dims, radius=np.pi,
                                          neuron_type=nengo.Direct(), label="output")
    # one ensemble array per output pattern
    # each ensemble array has the output dimensions
    # combine ensemble arrays to combine patterns
    # BONUS: figure out how to combine individual dimensions, once I know what they actually mean...
    # SUPER BONUS: use visual assement for the robot to be able to imitate a movement
    # ASIDE: which would be cool, because then maybe the robot could infer properties of objects from movement

    for n_i, nm in enumerate(pattern_file_names):
        name = nm[5:]
        print(name)
        n = make_dmp_net(function_list[n_i], osc[:2], output.input, name=name)
        dmp_net_list.append(n)

"""
    # probe the output
    p_out = nengo.Probe(output.output, synapse=0.003)

with nengo.Simulator(model) as sim:
    sim.run(4)

# un-normalise on export based off the original domain

plt.plot(sim.data[p_out])
plt.show()
"""