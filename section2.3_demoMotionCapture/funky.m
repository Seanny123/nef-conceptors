function [ x ] = funky( t )
% rise with three different linear slopes

t_scale = 0.5;

t_step = mod(t, 0.5);

if 0.1 * t_scale <  t_step < 0.3 * t_scale
    x = t_step - 0.1*t_scale;
elseif 0.3 * t_scale <= t_step < 0.6 * t_scale
    x = (0.3*t_scale)*(1+0.1) - t_step*0.1 - 0.1*t_scale;
elseif 0.6 * t_scale <= t_step
    x = (0.3*t_scale)*(1+0.1) - (0.6*t_scale)*(0.1+0.2) + 0.2*t_step - 0.1*t_scale;
else
    x = 0;
end

end

