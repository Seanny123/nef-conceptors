function d2 = mcresampleHJ(d, newfreq, method)
% Resamples motion capture data using interpolation.
%
% diff to original mcresample: field d.mus is also resampled
%
% syntax
% d2 = mcresample(d, newfreq, method);
%
% input parameters
% d: MoCap structure
% newfreq: new frame rate
% method: interpolation method (optional, default 'linear'; for other options, see help interp1)
%
% output
% d2: MoCap structure
%
% examples
% d2 = mcresample(d, 240);
% d2 = mcresample(d, 360, 'spline');
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland

d2 = [];
if isfield(d,'type') && (strcmp(d.type, 'MoCap data') || strcmp(d.type, 'norm data'))
    if nargin==2 
        method = 'linear'; 
    end
    d2 = d;
    t1 = (0:(size(d.data,1)-1))/d.freq;
    t2 = 0:(1/newfreq):t1(end);
    d2.data = interp1(t1, d.data, t2, method);
    d2.freq = newfreq;
    d2.nFrames = size(d2.data,1);
    if isfield(d,'mus') % addition HJ 
        d2.mus = (interp1(t1, d.mus', t2, method))';
    end
else
    disp('The first input argument should be a variable with MoCap data structure.');
end
