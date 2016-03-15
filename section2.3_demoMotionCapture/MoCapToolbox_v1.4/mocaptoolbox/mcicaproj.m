function [dp, p] = mcicaproj(d, pcs, ics)
% Performs an Independent Components analysis on MoCap, norm or segm data,  
% using the FastICA algorithm, and projects the data onto selected components.
%
% syntax
% [di, p] = mcicaproj(d, pcs, ics);
%
% input parameters
% d: MoCap, norm or segm data structure
% pcs: number of PCs entered into ICA
% ics: number of ICs estimated
%
% output
% di: vector of MoCap, norm or segm data structures
% p: ICA parameter structure containing the following fields:
%    icasig: independent components 
%    A: mixing matrix 
%    W: separation matrix 
%    meanx: mean vector of variables
%
% examples
% [di, p] = mcicaproj(d, 6, 3);
% 
% comments
% Uses the fastICA algorithm, implemented in the FastICA Package, which is available at http://www.cis.hut.fi/projects/ica/fastica/
% 
% see also
% mcpcaproj, mcsethares
% 
% © Part of the Motion Capture Toolbox, Copyright ©2008,
% University of Jyvaskyla, Finland

if isfield(d,'type') && (strcmp(d.type, 'MoCap data') || strcmp(d.type, 'norm data'))
    p = mcica(d, pcs, ics);
    for k=1:ics
        dp(k) = d;
        dp(k).data = repmat(p.meanx,d.nFrames,1)+(p.icasig(k,:)'*p.A(:,k)');
        dp(k).nFrames = size(dp(k).data,1);
    end


elseif isfield(d,'type') && strcmp(d.type, 'segm data')
    p = mcica(d, pcs, ics);
    for k=1:size(p.icasig,1) %% fixed 290509 PT
        dp(k) = d;
        tmp = repmat(p.meanx,d.nFrames,1)+(p.icasig(k,:)'*p.A(:,k)');
        dp(k).roottrans = tmp(:,1:3);
        for m=2:d.nMarkers
            dp(k).segm(m).eucl = tmp(:,3*m+(-2:0));
        end
        dp(k).nFrames = size(tmp,1);
    end
end
