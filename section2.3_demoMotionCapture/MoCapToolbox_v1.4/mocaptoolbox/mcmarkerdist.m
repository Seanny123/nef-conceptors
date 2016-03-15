function dist = mcmarkerdist(d, m1, m2)
% Calculates the frame-by-frame distance of a marker pair.
%
% syntax
% dist = mcmarkerdist(d, m1, m2);
%
% input parameters
% d: MoCap data structure
% m1, m2: marker numbers
%
% output
% dist: column vector
%
% examples
% dist = mcmarkerdist(d, 1, 5);
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland

if isfield(d,'type') && strcmp(d.type, 'MoCap data')
    if m1>d.nMarkers || m2>d.nMarkers || m1<1 || m2<1 
        disp('Marker numbers out of range.'); 
        return; 
    end
    c1 = 3*m1+(-2:0); c2 = 3*m2+(-2:0); 
    dist = sqrt(sum((d.data(:,c1)-d.data(:,c2)).^2,2));
else
    disp('The first input argument should be a variable with MoCap data structure.');
end

