function d2 = mcreorderdims(d, dims)
% Reorders the Euclidean dimensions in motion capture data.
%
% syntax
% d2 = mcreorderdims(d, dims);
%
% input parameters
% d: MoCap structure
% dims: vector containing the new order of dimensions
%
% output
% d2: MoCap structure
%
% examples
% d2 = mcreorderdims(d, [1 3 2]);
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland

d2 = [];
if isfield(d,'type') && strcmp(d.type, 'MoCap data')
    d2 = d;
    i2=[];
    for k=1:d.nMarkers
        i2=[i2 dims+3*k-3];
    end
    d2.data = d.data(:,i2);
else
    disp('The first input argument should be a variable with MoCap data structure.');
end
