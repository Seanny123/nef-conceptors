function d2 = mcs2j(d, par)
% Performs a segment-to-joint mapping.
%
% syntax
% d2 = mcs2j(d, par);
%
% input parameters
% d: segm data structure
% par: j2spar structure
%
% output
% d2: MoCap data structure
% 
% comments
% See the description of the j2spar structure.
% 
% see also
% mcinitj2spar, mcj2s
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland


d2=[];


if isfield(d,'type') && strcmp(d.type, 'segm data') && isfield(par,'type') && strcmp(par.type, 'j2spar')
    d2 = d;
    d2=rmfield(d2,'parent');
    d2=rmfield(d2,'roottrans');
    d2=rmfield(d2,'rootrot');
    d2=rmfield(d2,'segm');
    d2=rmfield(d2,'segmentName');
    d2.markerName = [];
    d2.type = 'MoCap data';
    d2.data = zeros(d2.nFrames, d2.nMarkers);
    d2.data(:, 3*par.rootMarker+(-2:0)) = d.roottrans;
    for k=1:d.nMarkers
        m=k;
        %tmp=d.roottrans;
        tmp=zeros(size(d.roottrans)); % added 050409 PT
        while(d.parent(m)>0)
            tmp = tmp + d.segm(m).eucl;
            m = d.parent(m);
        end
        d2.data(:,3*k+(-2:0)) = tmp;
    end
    d2 = mcrotate(d2, d.rootrot.az, [0 0 1], [0 0 0]); %% added 101208 PT
    d2.data = d2.data + repmat(d.roottrans,1,d.nMarkers); % added 050409 PT
else
    disp('The first input argument should be a variable with MoCap data structure.');
    disp('The second input argument should be a variable with j2spar data structure.');
end

return

