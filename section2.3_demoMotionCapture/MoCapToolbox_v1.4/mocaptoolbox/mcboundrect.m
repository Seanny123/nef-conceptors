function br = mcboundrect(d, mnum, w, hop)
% Calculates the bounding rectangle (the smallest rectangular area that contains the
% projection of the trajectory of each marker on the horizontal plane (i.e., floor).
%
% syntax
% br = mcboundrect(d);
% br = mcboundrect(d, mnum);
% br = mcboundrect(d, mnum, w, hop);
%
% input parameters
% d: MoCap data structure
% mnum: marker numbers
% w: length of analysis window (default: 4 sec)
% hop: overlap of analysis windows (default: 2 sec)
%
% output
% br: data matrix (windows x nMarkers) 
%
% examples
% br = mcboundrect(d);
% br = mcboundrect(d, [1 3 5]);
% br = mcboundrect(d, [1:d.nMarkers], 3, 1);
%
% comments
% If the function is called with the mocap data structure as the only input
% parameter, the calculation is performed for all markers with the default
% parameters. If the window and overlap length are to be changed, the
% markers have to be always specified (e.g., all markers by [1:d.nMarkers]).
%
% © Part of the Motion Capture Toolbox, Copyright ©2008,
% University of Jyvaskyla, Finland

if nargin==1
    w=4;
    hop=2;
end

if nargin>1
    d=mcgetmarker(d, mnum);
end

if nargin==2
    w=4;
    hop=2;
end

if nargin==3
    hop=2;
end

rect=[];

if isfield(d,'type') && (strcmp(d.type, 'MoCap data'))
    for k=1:d.nMarkers
        rtmp=[];
        for b=0:hop:(d.nFrames/120)-w
            ind1=1+120*b;
            ind2=min(size(d.data,1), ind1+120*w);
            tmp=d.data(ind1:ind2,k*3-2:k*3-1);%markers
            mintmp=min(tmp); 
            maxtmp=max(tmp);
            rtmp = [rtmp (maxtmp(1)-mintmp(1))*(maxtmp(2)-mintmp(2))/1000000];
        end
        rtmp=rtmp';
        rect = [rect rtmp];
    end
else 
    disp([10, 'This function only works with MoCap data structures.', 10]);
end

br = rect;
