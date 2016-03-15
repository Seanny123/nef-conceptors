function [amp, phase, h] = mchilbert(d)
% Calculates the Hilbert transform of data in a MoCap or norm structure.
% 
% syntax
% [amp, phase, h] = mchilbert(d);
% 
% input parameter
% d: MoCap or norm data structure
% 
% output
% amp: amplitude of analytic function derived from zero-mean signal
% phase: (unwrapped) phase of analytic function derived from zero-mean signal
% h: analytic function
%
% comments
% See help hilbert
% 
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland


if isfield(d,'type') && (strcmp(d.type, 'MoCap data') || strcmp(d.type, 'norm data'))
    amp=d; amp.type = 'norm data';
    phase=d; phase.type = 'norm data';
    h=d; h.type = 'norm data';
    m = mean(d.data); d.data = d.data - repmat(m,size(d.data,1),1);
    h.data = hilbert(d.data);
    amp.data = abs(h.data);
    phase.data = unwrap(angle(h.data));
    h.data = h.data + repmat(m,size(h.data,1),1);
else
    disp('The first input argument should be a variable with MoCap or norm data structure.');
end
