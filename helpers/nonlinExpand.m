function w = nonlinExpand(v, picks)
vv = kron(v,v);
vv = sign(vv) .* sqrt(abs(vv));
w = tanh(vv(picks));