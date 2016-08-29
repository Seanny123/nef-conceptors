import nengo

def make_good_osc(recurrent_tau, frequency, n_neurons, bump_func, neuron_type=nengo.LIFRate, seed=0):
    with nengo.Network() as test_net:
        osc = nengo.Network()
        osc.config[nengo.Ensemble].neuron_type = neuron_type
        osc.config[nengo.Ensemble].seed = seed
        nengo.networks.Oscillator(recurrent_tau, frequency, n_neurons, net=osc)
        bump = nengo.Node(bump_func)
        nengo.Connection(bump, osc.ensemble[0])

        p_osc = nengo.Probe(osc.ensemble, synapse=0.01)

    with nengo.Simulator(test_net) as test_sim:
        test_sim.run(10)

    zos = sim.trange()[np.isclose(np.zeros(dim.shape), dim, atol=0.001)]
    gain = np.mean(np.diff(zos)[2:]) * frequency

    good_osc = nengo.Network()
    good_osc.config[nengo.Ensemble].neuron_type = neuron_type
    good_osc.config[nengo.Ensemble].seed = seed
    nengo.networks.Oscillator(recurrent_tau, frequency*gain, n_neurons, net=good_osc)

    return good_osc


