from dmp_utils import *


def make_attrac_net(target_func, n_neurons=500, dd=None, seed=0, label=None):

    with nengo.Network(label=label) as ad_net:

        ad_net.input = nengo.Node(size_in=1)
        ad_net.output = nengo.Node(size_in=1)
        goal = nengo.Node([0])

        attractor = gen_point_attractor(ad_net, goal, n_neurons=n_neurons, seed=seed)
        nengo.Connection(attractor[0], ad_net.output, synapse=None)

        dest = target_func(np.linspace(-np.pi, np.pi, 100)).reshape((-1, 1))
        force_func = gen_forcing_functions(dest, num_samples=90)[0]

        if dd is not None:
            def conn_func(x, dec=dd, ff=force_func):
                return force_theta(x, dec, ff)
        else:
            def conn_func(x, ff=force_func):
                return force(x, ff)

        nengo.Connection(ad_net.input, attractor[1], synapse=None)

    return ad_net, conn_func
