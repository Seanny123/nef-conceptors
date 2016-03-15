function [dt, dn] = mcdecompose(d, order)
% Decomposes a kinematic variable into tangential and normal components.
%
% syntax
% [dt, dn] = mctangcomp(d, order);
%
% input parameters
% d: MoCap data structure containing either location or velocity data 
%   (timederorder = 0 or 1)
% order: time derivative order of the variable, must be at least 2 
%   (2 = acceleration, 3 = jerk, etc.)
%
% output
% dt: norm data structure containing the tangential components
% dn: norm data structure containing the normal components
%
% examples
% [dt, dn] = mcdecompose(d, 2); % acceleration
% [dt, dn] = mcdecompose(d, 3); % jerk 
% [dt, dn] = mcdecompose(d, 4); % jounce / snap
% [dt, dn] = mcdecompose(d, 5); % crackle
% [dt, dn] = mcdecompose(d, 6); % pop
% [dt, dn] = mcdecompose(d, 7); % you-name-it
%
% © Part of the Motion Capture Toolbox, Copyright ©2008,
% University of Jyvaskyla, Finland

if (order<2 || floor(order)~=order) disp('The second argument should be an integer greater than 1'); return; end

if isfield(d,'type') && strcmp(d.type, 'MoCap data')
   if d.timederOrder==0
       v = mctimeder(d, 1);
       q = mctimeder(d, order);
   elseif d.timederOrder==1
       v = d;
       q = mctimeder(d, order-1);
   else
       disp('The timederOrder field of the first argument shoud be 0 or 1'); return; 
   end
   vnorm = mcnorm(v);
   dt = vnorm;
   dn = vnorm;
   av =[];
   for k=1:3:(size(d.data,2))
       av = [av dot(v.data(:,(k:(k+2)))', q.data(:,(k:(k+2)))')']; % dot products per marker
   end
   dt.data = av ./ vnorm.data;
   anorm = mcnorm(q);
   dn.data = sqrt(anorm.data.^2 - dt.data.^2);
end

