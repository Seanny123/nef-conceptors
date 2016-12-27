function [ return_val ] = slow_triangle( t )
% rise with three different linear slopes
persistent x;

mag_scale = 0.002;
t_scale = 0.5;

t_step = mod(t, t_scale);

if t_step == 0
    x = 0.1;
end

if t_step < 0.5*t_scale
    x = x + (1 * mag_scale);
elseif 0.5*t_scale <= t_step
    x = x - (1 * mag_scale);
else
    x = 0.1;
end

return_val = x;

end

