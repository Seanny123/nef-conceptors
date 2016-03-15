function p = mcanimate_testDiscard(d, p, proj)
% Creates png frames for an animation of motion capture data, and saves them to a folder.
%
% syntax
% p = mcanimate(d);
% p = mcanimate(d, p);
% p = mcanimate(d, p, proj);
% 
% input parameters
% d: MoCap data structure
% p: animpar structure
% proj: projection used: 0 = orthographic (default), 1 = perspective
%
% output
% p: animpar structure used for plotting the frames
% 
% examples
% mcanimate(d, par, 1);
%
% comments
% If the animpar structure is not given as input argument, the function
% creates it by calling the function mcinitanimpar and setting the .limits field 
% of the animpar structure automatically so that all the markers fit into all frames.
% If the par.pers field (perspective projection) is not given, it is created internally for backwards
% compatibility. For explanation of the par.pers field, see help mcinitanimpar
%
% see also
% mcplotframe, mcinitanimpar
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland


if nargin<3
    proj=0;
end
if nargin<2
    p = mcinitanimpar;
end

% resample for p.fps fps movie
d2 = mcresample(d, p.fps);

p.animate = 1;

%for compatibility
if ~isfield(p, 'pers')
    p.pers.c=[0 -3000 0];
    p.pers.th=[0 0 0];
    p.pers.e=[0 -2000 0];
end

p = mcplotframe_testDiscard(d2,1:d2.nFrames,p, proj); % output parameter added 240608

