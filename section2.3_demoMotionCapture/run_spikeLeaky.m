addpath('../helpers');

leakRates = [0 0.001 0.01];
leakNum = length(leakRates);

tMax = 0.5;
tSteps = 0:0.001:tMax;
tLen = length(tSteps);

p3 = zeros(size(tSteps'));
tmp_count = 1;
for t = tSteps
   p3(tmp_count) = triangle(t);
   tmp_count = tmp_count + 1;
end
d = fdesign.lowpass('N,Fc', 30, 5, 1/0.001);
Hd = design(d);
f_p3 = filter(Hd,p3);

sig = zeros(tLen - 50);
testLen = 2000;
results = zeros(leakNum, testLen);

for l_i = 1:leakNum
    [sig, results(l_i, :)] = varLeaky(f_p3, leakRates(l_i));
end

legendCell = cell(1, leakNum+1);
hold on;
tSteps = 0:0.001:((testLen+50)*0.001);
refSig = zeros(size(tSteps'));
tmp_count = 1;
for t = tSteps
   refSig(tmp_count) = triangle(t);
   tmp_count = tmp_count + 1;
end;

f_ref = filter(Hd, refSig);

plot(f_ref(50:end), 'k');
legendCell{1} = 'reference';
for l_i = 1:leakNum
    plot(results(l_i, :));
    legendCell{1+l_i} = sprintf('LR=%0.2f', leakRates(l_i));
end

legend(legendCell);
hold off;
