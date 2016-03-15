function dj = nnRaw2jHJ_absH(nnRawData, hmean, segLengths, j2spar, freq)
% transform nnRawData back to mocap structure. Use raw Data where root
% height is given in absolute values

% extract components from nnRawData
hraw = nnRawData(:,1);
forwardAngDiffsraw = nnRawData(:,2);
xyDiffRelativeraw = nnRawData(:,3:4);
segUnitVecsraw = nnRawData(:,5:end);

nSegments = (size(nnRawData,2) - 4) / 3;
N = size(nnRawData,1); % nr of frames

% initialize joint cartesian data
nFrames = size(nnRawData, 1);
jdata = zeros(nFrames, 3 * (nSegments + 1));

%%% compute root trajectory

% fill root height
jdata(:,3) = hraw;  % assuming the root index is 1
% integrate body normal vec xy angle
forwardAngAbs = cumsum(forwardAngDiffsraw);
%%% fill root xy absolute motion
% rotate xyDiffRelativeraw by forwardAngAbs
[theta, r] = cart2pol(xyDiffRelativeraw(:,1), xyDiffRelativeraw(:,2));
[xyDiffAbs1, xyDiffAbs2] = pol2cart(theta + forwardAngAbs, r);
xyAbs = [cumsum(xyDiffAbs1), cumsum(xyDiffAbs2)];
jdata(:,1:2) = xyAbs;

%%% fill remaining joint trajectories
% rotate segUnitVecsraw by forwardAngAbs
euclidsRot = cell(1, nSegments);
for i = 1:nSegments
    euclidsRot{i} = segLengths(i) * ...
        rotatexyHJ(segUnitVecsraw(:,(i-1)*3+1:i*3),forwardAngAbs);
end
% sum along kinematic chain
parents = j2spar.parent;
for k = 1:nSegments
    m = k + 1;
    tmp = zeros(N, 3);
     while(parents(m) > 0)
        tmp = tmp + euclidsRot{m - 1};
        m = parents(m);
     end
    tmp = tmp + jdata(:,1:3);   
    jdata(:,3*k+(1:3)) = tmp;
end

% create MoCap structure
dj = mcinitstruct( 'MoCap data', jdata, freq);
dj.nFrames = N;




