function i = randpick(p)
% picks a index from a probability vector p according to the probs in p
if abs(sum(p) - 1) > 1e-14 
    error('randpick needs prob vector argument');
end
pcs = cumsum(p);
rd = rand;
i = sum(pcs <= rd) + 1;
