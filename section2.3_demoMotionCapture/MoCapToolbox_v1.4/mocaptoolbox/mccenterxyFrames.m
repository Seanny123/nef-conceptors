function d2 = mccenterxyFrames(d)
% Translates motion capture data to have a centroid of [0 0] in x,y across markers per frame.
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

d2=d;

L = size(d.data,1); % runlength
nJ = size(d.data,2) / 3; % nr of joints
for n = 1:L

    x = mean(d.data(n,1:3:end));
    y = mean(d.data(n,2:3:end));
    d2.data(n,1:3:end) = d.data(n,1:3:end) - x * ones(1,nJ);
    d2.data(n,2:3:end) = d.data(n,2:3:end) - y * ones(1,nJ);
    
end
