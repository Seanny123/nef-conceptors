function [nnRawData, hmean, segLengths] = j2nnRawHJ(d, par)
% transform mocap structure d to nn raw data
%
% the mocap data matrix must have the root in first three columns, the left
% hip in cols 16-18, the right hip in cols 4-6
%
% Inputs: d: Mocap structure
%         par: j2spar data structure
%
% Outputs:
% nnRawData: matrix size nFrames x (4 + 3*nSegments) = nFrames x 81 which
%            is organized as follows:
%            col1: hDiffs
%            col2: forwardAngDiffs
%            cols 3&4: xyDiffsRelative
%            remaning cols: segUnitVecs
% hmean: mean value of root height above ground
% segLenghts: vector size 1 x nSegments, giving the mean (over frames) of
%            segment lengths

N = size(d.data, 1); % length of recording (Nr of frames)
Nj = size(d.data,2) / 3; % nr of joints

roottrans = d.data(:,1:3);

% shift joint data to zero root
jDataZeroRoot = d.data - repmat(roottrans, 1, Nj);

lHipZeroRoot = jDataZeroRoot(:,16:18); % lhip is marker # 6
rHipZeroRoot = jDataZeroRoot(:,4:6); % lhip is marker # 2

hipVecZeroRoot = lHipZeroRoot - rHipZeroRoot; % vec from rhip to lhip
hipVecZeroRootxy = hipVecZeroRoot(:,1:2);
hipVecZeroRootxyNorms = sqrt(sum(hipVecZeroRootxy.^2,2));
hipVecZeroRootxyUnit = diag(1 ./ hipVecZeroRootxyNorms) * hipVecZeroRootxy;
% rotate by -pi/2
[theta r] = cart2pol(hipVecZeroRootxyUnit(:,1),...
    hipVecZeroRootxyUnit(:,2));
rootrotAbs = theta - pi/2; % xy rotation angle of body forward vector
[X, Y] = pol2cart(rootrotAbs, r);
forwardxyUnit = [X, Y]; % unit xy vector of forward direction

%%% collect root motion data
%  height differentials hDiffs
hraw = roottrans(:,3); % height of root above ground
hDiffs = zeros(N,1);
hDiffs(1,1) = hraw(2) - hraw(1);
hDiffs(2:end,1) = hraw(2:end,1) - hraw(1:end-1, 1); %COLLECT

hmean = mean(hraw);

% angular diffs of forwardxyUnit
forwardAngDiffs = zeros(N,1);
rootrotAbsUnwrapped = unwrap(rootrotAbs);
forwardAngDiffs(2:end, :) = rootrotAbsUnwrapped(2:end,:) - ...
    rootrotAbsUnwrapped(1:end-1,:);
forwardAngDiffs(1, :) = forwardAngDiffs(2, :); %COLLECT

% xy Diffs relative to forwardxyUnit
xyraw = roottrans(:,1:2);
xyDiffsAbs = zeros(N,2);
xyDiffsAbs(1,:) = xyraw(2,:) - xyraw(1,:);
xyDiffsAbs(2:end,:) = xyraw(2:end,:) - xyraw(1:end-1,:); % in absolute xy
[xyDiffsTheta xyDiffsr] = cart2pol(xyDiffsAbs(:,1), xyDiffsAbs(:,2));
[xyDiffsRelativex, xyDiffsRelativey]  = ...
    pol2cart(xyDiffsTheta - rootrotAbs, xyDiffsr);
xyDiffsRelative = [xyDiffsRelativex, xyDiffsRelativey];  % COLLECT

%%% collect angular data of non-root joints (assuming root = joint #1)
% rotate root=zero data by -rootrotAbs
jDataZeroRootRot = jDataZeroRoot;
for j = 2:Nj
    jDataZeroRootRot(:,(j-1)*3+1:j*3) = ...
        rotatexyHJ(jDataZeroRoot(:,(j-1)*3+1:j*3), -rootrotAbs);
end
% collect segment unit vectors
segUnitVecs = zeros(N, (Nj - 1)*3); % COLLECT
segLengths = zeros(1, Nj-1); % COLLECT
% Euclidean parameters
for k = 1:length(par.parent)
    if par.parent(k)>0
        ind1 = 3*k + (-2:0);
        ind2 = 3*par.parent(k) + (-2:0);
        segVec = jDataZeroRootRot(:,ind1) - jDataZeroRootRot(:,ind2);
        segVecNorms = sqrt(sum(segVec.^2,2));
        segLengths(1,k-1) = mean(segVecNorms);
        segVecUnit = diag(1 ./ segVecNorms) * segVec;
        segUnitVecs(:, 3*(k-1) + (-2:0)) = segVecUnit;
    end
end
nnRawData = [hDiffs forwardAngDiffs xyDiffsRelative segUnitVecs];
