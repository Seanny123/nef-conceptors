function newmus  = smoothmus( mus )
dim = size(mus,1); L = size(mus,2);
newmus = mus;
for n = 1:L
    for d = 1:dim
        if mus(d,n) >= 0.5 && mus(d,n) ~= 1.0
            newmus(d,n) = 1 - 2 * (1 - mus(d,n))^2;
        elseif mus(d,n) < 0.5 && mus(d,n) ~= 0
            newmus(d,n) = 2 * (mus(d,n))^2;
        end
    end
end

end

