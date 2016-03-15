function d2 = mcnorm(d, comps)
% Calculates the norms of kinematic vectors.
%
% syntax
% n = mcnorm(d);
% n = mcnorm(d, comps);
%
% input parameters
% d: MoCap data structure
% comps: components included in the calculation (optional, default = 1:3)
%
% output
% n: norm data structure
%
% examples
% n = mcnorm(d);
% n = mcnorm(d, 1:2); % calculates norm of horizontal projection
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland

if nargin==1 
    comps=1:3; 
end
d2=[];
if isfield(d,'type') && strcmp(d.type, 'MoCap data')
    d2 = d;
    d2.type = 'norm data';
    d2.data=[];
    for k=1:3:size(d.data,2)
        d2.data = [d2.data sqrt(sum(d.data(:,k+comps-1).^2,2))];
    end
else
    d2=[];
    for k=1:3:size(d,2)
        d2 = [d2 sqrt(sum(d(:,k+comps-1).^2,2))];
    end

end
