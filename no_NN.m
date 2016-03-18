%%%% Main script for motion capture demo
%%%% H. Jaeger, May 23, 2015
%%%%
%%%% External resources included:
%%%% Motion capture raw data from http://mocap.cs.cmu.edu/
%%%% Motion capture analysis and visualization Matlab toolbox from
%%%% https://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/mocaptoolbox
%%%%
%%%% Note: users running later versions of Matlab reported that some
%%%% conceptor demos wouldn't succeed. The cause seems to be that my old
%%%% Matlab version uses a different (and better) implementation of the inv
%%%% function than later Matlab versions. Replacing all calls to "inv" by
%%%% calls to "pinv" so far has always resolved the problems. 

addpath('./MoCapToolbox_v1.4/mocaptoolbox');
addpath('./nnData');
addpath('./helpersMocap');
addpath('../helpers');
addpath('./makeNNData'); % comment out if you don't want to 
                         % preprocess raw mocap data

set(0,'DefaultFigureWindowStyle','docked');

%%% Experiment control
randstate = 1; % random seed for creating network weights
newNets = 1; % whether new raw nets are created
newSystemScalings = 1; % whether raw nets are freshly re-scaled
reloadData = 1; % whether normalized neural-network input data are
% re-loaded
relearn = 1; % whether model is re-trained
makeAllDataFlag = 0; % whether normalized neural-network input data are
% re-computed from downloaded mocap source data
% (expensive)
showStates = 0; % whether to plot some reservoir states obtained during testing
showSingVals = 0; % whether to plot singular value spectra of conceptors
showOutTraces = 0; % whether to plot outputs of some of the 61 pattern dimensions
                   % collected for the entire demo run
showVideo = 1; % whether video is generated and shown (Matlab-internal, slow,
% does not create an externally runnable video in one of the
% standard video formats)

%%% a few constants
nP = 15; % nr of patterns; if you experiment with additional patterns
% (or subsets of the patterns used here), adapt
pattDim = 61; % dimension of patterns, don't change from value 61!


%%% Setting system params
Netsize = 600; % network size
NetSR = 1; % spectral radius
NetinpScaling = .8 * ones(pattDim, 1); % scaling of input weights
NetinpScaling([5 17],1) = 0 * ones(2,1); % these two were found to "function"
% primarily as noise and are
% suppressed
BiasScaling = .8; % scaling of bias vector
LR = .6; % leaking rate for leaky-integration neuron model, ranges in [0 1]
% with a value of 1 meaning no integration -- the smaller, the
% more smoothing

washoutLength = 50; % how many initial updates per pattern are discarded
% in loading and Wout computation

%%% Ridge regression regularization. Set to 0 because data were so
%%% variable that no regularization was found necessary
TychonovAlpha = 0; % regularizer for  W training
TychonovAlphaReadout = 0; % regularizer for  Wout training

%%% Conceptor learning and testing
alphas = 10 * ones(1,nP); % apertures
alphas(2) = 3;
alphas(13) = 20;
CtestLength = 1000;

%%% numbering of motion patterns:
% 1 ExaggeratedStride 2 SlowWalk 3 Walk 4 RunJog 5 CartWheel  6 Waltz
% 7 Crawl  8 Standup  9 Getdown 10 Sitting  11 GetSeated  12 StandupFromStool
% 13 Box1  14 Box2 15 Box3

%%% setting sequence order of patterns to be displayed
pattOrder = [10 12 2 1 4 1 6 3 9 7 8 3 13 15 14 2  5 3 2 11 ];
%%% setting durations of each pattern episode (note: 120 frames = 1 sec)
pattDurations = [ 150 260 200 200 250 130 630 200 120 400 100 100  ...
    250 400 300 150 670 100 150 300 ];
%%% setting durations for morphing transitions between two subsequent patterns
pattTransitions = 120 * ones(1, length(pattOrder)-1);

plotPicks = [3 4 25 26  28]; % which input dimensions are plotted

%%% additive state noise level inserted during testing and video generation
stateNL = 0.;

if makeAllDataFlag
    makeAllData;
end


if reloadData
    load nnRawExaStride; load nnRawSlowWalk;  load nnRawWalk;
    load nnRawRunJog; load nnRawCartWheel; load nnRawWaltz;
    load nnRawCrawl; load nnRawStandup; load nnRawGetdown;
    load nnRawSitting; load nnRawGetSeated; load nnRawStandupFromStool;
    load nnRawBox1; load nnRawBox2;  load nnRawBox3;
    
    
    p1 = nnRawDataExaStride; p2 = nnRawDataSlowWalk;
    p3 = nnRawDataWalk; p4 = nnRawDataRunJog;
    p5 = nnRawDataCartWheel; p6 = nnRawDataWaltz;
    p7 = nnRawDataCrawl; p8 = nnRawDataStandup;
    p9 = nnRawDataGetdown; p10 = nnRawDataSitting;
    p11 = nnRawDataGetSeated; p12 = nnRawDataStandupFromStool;
    p13 = nnRawDataBox1;  p14 = nnRawDataBox2; p15 = nnRawDataBox3;
    
    % set segment lenghts of stick figure to the ones found in the
    % GetSeated data
    segmentlengths = segLengthsGetSeated;
    
    % set default height of center of gravity to mean of some of the traces
    % (needed to place visualized stick guy in a reasonably looking height
    % above ground)
    hmean = mean([hmeanExaStride, hmeanSlowWalk, ...
        hmeanWalk]);
    
    % pattern durations
    pl1 = size(p1,1); pl2 = size(p2,1); pl3 = size(p3,1);
    pl4 = size(p4,1); pl5 = size(p5,1); pl6 = size(p6,1);
    pl7 = size(p7,1); pl8 = size(p8,1); pl9 = size(p9,1);
    pl10 = size(p10,1); pl11 = size(p11,1); pl12 = size(p12,1);
    pl13 = size(p13,1); pl14 = size(p14,1); pl15 = size(p15,1);
    
    pattlengths = [pl1 pl2 pl3 pl4 pl5 pl6  pl7  pl8 ...
        pl9 pl10 pl11 pl12 pl13 pl14 pl15];
    
    % normalize all raw data simultaneously to range -1 to +1
    allData = [p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11; p12; ...
        p13; p14; p15];
    [nnNormAllData, scalings, shifts] = ...
        normalizeDataMinus1Plus1(allData);
    
    % back-distribute concatenated traces over the individual patterns
    patts = cell(1,nP);
    startInd = 1;
    for i = 1:nP
        patts{i} = ...
            nnNormAllData(startInd:startInd + pattlengths(i)-1,:);
        startInd = startInd + pattlengths(i);
    end
end

%%% Initializations

randn('state', randstate);
rand('twister', randstate);

%%% create video
if showVideo
    j2spar.type = 'j2spar';
    j2spar.rootMarker = 1;
    j2spar.frontalPlane = [6 10 2 ];  % be careful about order
    j2spar.parent = [0 1 2 3 4 1 6 7 8 1 10 11 11 13 14 15 11 17 18 19];
    j2spar.segmentName = cell(1,19);
    
    freq = 120;
    
    djrecovered = ...
        nnRaw2jHJ_absH(p1, hmean, segmentlengths, ...
        j2spar, freq);
    
    japarVideo.showmnum = 0;
    japarVideo.animate = 1;
    japarVideo.colors = 'wkkkk';
    japarVideo.trl = 0;
    japarVideo.numbers = [];
    japarVideo.cwidth = 5;
    japarVideo.twidth = 1;
    japarVideo.az = 20;
    japarVideo.el = 20;
    japarVideo.fps = 30;
    japarVideo.limits = [];
    japarVideo.scrsize = [400 350];
    %japarVideo.scrsize = [400 300];
    japarVideo.showfnum = 0;
    japarVideo.conn2 = [];
    japarVideo.conn = [11 12; 13 11; 11 17; 13 14; 14 15; 15 16;...
        17 18; 18 19; 19 20; 10 11; 1 10; 1 2; 1 6; 2 3; 3 4; 4 5;...
        6 7; 7 8; 8 9];
    japarVideo.msize = 12;
    japarVideo.fontsize = 12;
    japarVideo.folder = '0'; % set to string '0' if no save is wanted, else
    % give a name for a folder where frame png
    % files are to be stored
    japarVideo.center = 1;
    japarVideo.pers.c = [0 -10000 0];
    japarVideo.pers.th = [0 0 0];
    japarVideo.pers.e = [0 -6000 0];
    japarVideo.nGridlines = 5;
    japarVideo.botWidth = 5000;
    
    japarVideo.fignr = 11;
    djrecoveredResampled = mcresampleHJ(djrecovered, japarVideo.fps);
    xRootShifts = djrecoveredResampled.data(2:end,1) - ...
        djrecoveredResampled.data(1:end-1,1);
    yRootShifts = djrecoveredResampled.data(2:end,2) - ...
        djrecoveredResampled.data(1:end-1,2);
    
    djrecoveredCentered = mccenterxyFrames(djrecovered);
    djrecoveredCenteredResampled = ...
        mcresampleHJ(djrecoveredCentered, japarVideo.fps);
    
    %%s part of the mocap toolbox
    mcanimateHJ(djrecoveredCenteredResampled, japarVideo, 1,...
        xRootShifts, yRootShifts);
end
