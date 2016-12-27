tMax = 1.5;
tSteps = 0:0.001:tMax;

p3 = zeros(size(tSteps'));

tmp_x = 0;
tmp_count = 1;
for t = tSteps
   p3(tmp_count) = triangle(t);
   tmp_count = tmp_count + 1;
end

d = fdesign.lowpass('N,Fc', 30, 5, 1/0.001);
Hd = design(d);
f_p3 = filter(Hd,p3);
hold on
plot(p3)
plot(f_p3);
hold off
