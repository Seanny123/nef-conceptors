addpath('../helpers');

leakRates = 0.1:0.1:0.3;
leakNum = length(leakRates);
sinPer = (2 * pi * 10);

tMax = 0.5;
tSteps = 0:0.001:tMax;
p1 = sin(tSteps*sinPer)';
noise = rand(size(p1));
p1 = p1 + noise/3;
tLen = length(tSteps);

sig = zeros(tLen - 50);
results = zeros(leakNum, tLen - 50);

for l_i = 1:leakNum
    [sig, results(l_i, :)] = noRec(p1, leakRates(l_i));
end

legendCell = cell(1, leakNum+1);
hold on;
plot(sig);
legendCell{1} = 'reference';
for l_i = 1:leakNum
    plot(results(l_i, :));
    legendCell{1+l_i} = sprintf('LR=%0.1f', leakRates(l_i));
end

legend(legendCell);
hold off;
