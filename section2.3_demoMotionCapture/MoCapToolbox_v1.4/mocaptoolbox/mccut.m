function [d1, d2] = mccut(d1, d2)
% Cuts two MoCap structures to the length of the shorter one.
%
% syntax
% [d1, d2] = mccut(d1, d2);
%
% input parameters
% d1, d2: MoCap or norm structures
%
% output
% d1, d2: MoCap or norm structures, one shortened and one original (both with same number of frames)
%
% © Part of the Motion Capture Toolbox, Copyright ©2008,
% University of Jyvaskyla, Finland

if d1.nFrames<d2.nFrames
    disp('d2 has been shortened.');
elseif d1.nFrames>d2.nFrames
    disp('d1 has been shortened.');
end

if d1.nFrames ~= d2.nFrames
    N = min(d1.nFrames, d2.nFrames);
    d1.data = d1.data(1:N,:);
    d2.data = d2.data(1:N,:);
    d1.nFrames = N;
    d2.nFrames = N;
end
