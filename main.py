import nengo
import numpy as np

T = 1.0
dt = 0.001

t = np.arange(1000)*dt

pattern1 = 2 * t -1
pattern2 = np.sin(t*2*np.pi)
pattern3 = np.where(t<0.5, 2*t, 2-2*t)
pattern4 = -pattern3

def make_pattern():
    sig = nengo.processes.WhiteSignal(T, high=2)
    # this is returning a flat signal and I don't know why
    return sig.run(1000)

np.random.seed(3)
#patterns = [make_pattern()[:,0] for i in range(4)]
patterns = [pattern1, pattern2, pattern3, pattern4]
# maps from input value (in this case, theta) to output value
patterns = np.array(patterns)

model = nengo.Network()
#model.config[nengo.Ensemble].neuron_type = nengo.LIFRate()
#model.config[nengo.Connection].solver = nengo.solvers.LstsqL2(reg=0.001)
tau = 0.1
omega = 2*np.pi
with model:
    system = nengo.Ensemble(n_neurons=1000, dimensions=3, radius=1.4, label='cycle')

    def cycle(x):
        """makes a speed controlled oscillator"""
        a = 1
        b = 2 * np.pi * x[2]
        r = np.sqrt(x[0]**2 + x[1]**2)
        theta = np.arctan2(x[1], x[0])
        dr = 10*(-r**3 + a*r)
        dtheta = b

        dx = dr*np.cos(theta) - r*np.sin(theta)*dtheta
        dy = dr*np.sin(theta) + r*np.cos(theta)*dtheta

        return [x[0] + tau*dx, x[1] + tau*dy]

    nengo.Connection(system, system[:2], synapse=tau, function=cycle)

    rate = nengo.Node([1])
    nengo.Connection(rate, system[2])

    readouts = nengo.networks.EnsembleArray(n_neurons=30, n_ensembles=len(patterns))
    for i in range(len(patterns)):
        def read_func(x, index=i):
            """based off angle, read out pattern"""
            theta_p = np.arctan2(x[1], x[0])
            r_p = np.sqrt(x[0]**2 + x[1]**2)
            if theta_p < 0: 
                theta_p += 2 * np.pi
            index_p = int(len(patterns[0]) * theta_p / (2 * np.pi))

            return patterns[index][index_p] * r_p
        nengo.Connection(system, readouts.ensembles[i], function=read_func)


    select = nengo.Node([2.5] * len(patterns))

    for i in range(len(patterns)):
        nengo.Connection(select[i], readouts.ensembles[i].neurons, transform=[[-1]]*30)

    readout = nengo.Ensemble(n_neurons=30, dimensions=1)
    nengo.Connection(readouts.output, readout, transform=[[1]*len(patterns)])