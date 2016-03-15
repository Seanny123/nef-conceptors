function d2 = mccenter(d)
% Translates motion capture data to have a centroid of [0 0 0] across markers and over time.
%
% syntax
% d2 = mccenter(d);
%
% input parameters
% d: MoCap structure or data matrix
%
% output
% d2: MoCap structure or data matrix
%
% comments
% Missing data (NaN's) is ignored when calculating the centroid.
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland

d2=[];
if isfield(d,'type') && strcmp(d.type, 'MoCap data')
    x = mean(mcmean(d.data(:,1:3:end)));
    y = mean(mcmean(d.data(:,2:3:end)));
    z = mean(mcmean(d.data(:,3:3:end)));
    d2 = mctranslate(d, [-x -y -z]);
else
    x = mean(mcmean(d(:,1:3:end)));
    y = mean(mcmean(d(:,2:3:end)));
    z = mean(mcmean(d(:,3:3:end)));
    d2 = mctranslate(d, [-x -y -z]);
end

