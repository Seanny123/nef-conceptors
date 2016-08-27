function [ return_val ] = so_funky( t )
% rise with three different linear slopes
persistent x;

if t == 0
    x = 0;
end

mag_scale = 0.002;
t_scale = 0.5;

t_step = mod(t, 0.5);

if 0.1*t_scale <  t_step && t_step < 0.3*t_scale
    x = x + (1 * mag_scale);
elseif 0.3*t_scale <= t_step && t_step < 0.6*t_scale
    x = x - (0.1 * mag_scale);
elseif 0.6*t_scale <= t_step
    x = x + (1.3 * mag_scale);
else
    x = 0;
end

return_val = x;

end

