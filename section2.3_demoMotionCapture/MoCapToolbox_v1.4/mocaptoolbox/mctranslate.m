function d2 = mctranslate(d, X)
% Translates motion-capture data by a vector.
%
% syntax
% d2 = mctranslate(d, X);
%
% input parameters
% d: MoCap structure or data matrix
% X: translation vector
%
% output
% d2: MoCap structure or data matrix
%
% examples
% d2 = mctranslate(d, [0 1000 0]);
%
% © Part of the Motion Capture Toolbox, Copyright ©2008, 
% University of Jyvaskyla, Finland

d2=[];
if isfield(d,'type') && strcmp(d.type, 'MoCap data')
    d2 = d;
    d2.data(:,1:3:end) = d2.data(:,1:3:end) + X(1);
    d2.data(:,2:3:end) = d2.data(:,2:3:end) + X(2);
    d2.data(:,3:3:end) = d2.data(:,3:3:end) + X(3);
else
    d2 = d;
    d2(:,1:3:end) = d(:,1:3:end) + X(1);
    d2(:,2:3:end) = d(:,2:3:end) + X(2);
    d2(:,3:3:end) = d(:,3:3:end) + X(3);
end

