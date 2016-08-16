import nengo

SEED = 0

# for each sig


# get the rate data
with model as nengo.Network:
    in_sig = nengo.Node()
    reserv = nengo.Ensemble(300, 1, neuron_type=nengo.LIFRate(), seed=SEED)

# do SVD on the neuron data

# use it for Conceptor testing
with model as nengo.Network:
    reserv = nengo.Ensemble(300, 1, neuron_type=nengo.LIFRate(), seed=SEED)
    

