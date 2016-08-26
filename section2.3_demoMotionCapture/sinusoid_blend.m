addpath('../helpers');

set(0,'DefaultFigureWindowStyle','docked');

%%% Experiment control
randstate = 1; % random seed for creating network weights
showStates = 1; % whether to plot some reservoir states obtained during testing
showSingVals = 1; % whether to plot singular value spectra of conceptors

%%% a few constants
nP = 2; % nr of patterns; if you experiment with additional patterns
% (or subsets of the patterns used here), adapt
pattDim = 1; % dimension of patterns, don't change from value 61!


%%% Setting system params
Netsize = 600; % network size
NetSR = 1; % spectral radius
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
CtestLength = 1000;

%%% setting sequence order of patterns to be displayed
pattOrder = [1 2];
%%% setting durations of each pattern episode
tMax = 0.5;
tSteps = 0:0.001:tMax;
tLen = size(tSteps, 2);
pattDurations = [tLen, tLen];
%%% setting durations for morphing transitions between two subsequent patterns
pattTransitions = 100 * ones(1, length(pattOrder)-1);

plotPicks = [1 1 1 1 1]; % which input dimensions are plotted


%%% Initialise the patterns
sinPer = (2 * pi * 10) / tMax;
cosPer = (2 * pi * 20) / tMax;
p1 = sin(tSteps*sinPer)';
p2 = 0.5*cos(tSteps*cosPer)';

p3 = zeros(tLen);

tmp_x = 0;
tmp_count = 1;
for t = tSteps
   p3(tmp_count) = so_funky(t);
   tmp_count = tmp_count + 1;
end

% pattern durations
pl1 = size(p1,1); pl2 = size(p2,1);
pattlengths = [pl1 pl2];
allData = [p1; p2;];

% back-distribute concatenated traces over the individual patterns
patts = cell(1,nP);
startInd = 1;
for i = 1:nP
    patts{i} = allData(startInd:startInd + pattlengths(i)-1,:);
    startInd = startInd + pattlengths(i);
end

%%% Initializations

randn('state', randstate);
rand('twister', randstate);


% Create raw weights
if Netsize <= 20
    Netconnectivity = 1;
else
    Netconnectivity = 10/Netsize;
end
WstarRaw = generate_internal_weights(Netsize, Netconnectivity);
WinRaw = 2*(rand(Netsize, pattDim) - 0.5);
WbiasRaw = 2*(rand(Netsize, 1) - 0.5);

% Scale raw weights and initialize weights
Wstar = NetSR * WstarRaw;
Win =  WinRaw;
Wbias = BiasScaling * WbiasRaw;

% harvest data from network externally driven by patterns
totalDataLength = sum(pattlengths);
totalLearnLength = totalDataLength - nP * washoutLength;

allTrainxArgs = zeros(Netsize, 0);
allTrainOldxArgs = zeros(Netsize, 0);
allTrainWtargets = zeros(Netsize, 0);
allTrainOuts = zeros(pattDim, 0);
Wtargets = zeros(Netsize,0);
patternRs = cell(1,nP);
startXs = zeros(Netsize, nP);
filename = '';
% collect data from driving native reservoir with different drivers
for p = 1:nP
    patt = patts{p}; % current pattern
    learnLength = pattlengths(p) - washoutLength;
    xCollector = zeros(Netsize, learnLength );
    xOldCollector = zeros(Netsize, learnLength );
    WTargetCollector = zeros(Netsize, learnLength);
    pCollector = zeros(pattDim, learnLength );
    x = zeros(Netsize, 1);
    for n = 1:(washoutLength + learnLength)
        u = patt(n,:)'; % pattern input
        xOld = x;
        Wtarget = Wstar * xOld + Win * u;
        x = (1-LR)*xOld + LR * tanh(Wtarget + Wbias);
        if n == washoutLength
            startXs(:,p) = x;
        end
        if n > washoutLength
            xCollector(:, n - washoutLength ) = x;
            xOldCollector(:, n - washoutLength ) = xOld;
            WTargetCollector(:, n - washoutLength ) = Wtarget;
            pCollector(:, n - washoutLength) = u;

        end
        uOld = u;
    end
    %filename = sprintf('patt%d.mat', p);
    %save(filename, 'xCollector');
    
    patternRs{p} = xCollector * xCollector' / learnLength;

    allTrainxArgs = [allTrainxArgs, xCollector];
    allTrainOldxArgs = [allTrainOldxArgs, xOldCollector];
    allTrainOuts = [allTrainOuts, pCollector];
    allTrainWtargets = [allTrainWtargets, WTargetCollector];
end

%save('corrs', 'patternRs')

%%% compute pattern readout
Wout = (pinv(allTrainxArgs * allTrainxArgs' + ...
    TychonovAlphaReadout * eye(Netsize)) ...
    * allTrainxArgs * allTrainOuts')';
% training error
outsRecovered = Wout*allTrainxArgs;
NRMSE_readout = nrmse(outsRecovered, allTrainOuts);
disp(sprintf('mean NRMSE readout: %g   mean abs Wout: %g',...
    mean(NRMSE_readout(not(isnan(NRMSE_readout)),1)),...
    mean(mean(abs(Wout)))));

outsTrain = cell(1,nP);
wo = washoutLength;
startInd = 1;
for i = 1:nP
    outsTrain{i} = outsRecovered(:,...
        startInd:startInd + pattlengths(i) - wo -1)';
    startInd = startInd + pattlengths(i) - wo;
end


%%% compute W
W = (pinv(allTrainOldxArgs * allTrainOldxArgs' + ...
    TychonovAlpha * eye(Netsize)) * allTrainOldxArgs * allTrainWtargets')';
% training errors per neuron
NRMSE_W = nrmse(W*allTrainOldxArgs, allTrainWtargets);
disp(sprintf('mean NRMSE W: %g   mean abs W: %g ', ...
    mean(NRMSE_W), mean(mean(abs(W)))));

%%% compute conceptors
Cs = cell(4, nP);
for p = 1:nP
    [U S V] = svd(patternRs{p});
    Snew = (S * pinv(S + alphas(p)^(-2) * eye(Netsize)));
    C = U * Snew * U';
    Cs{1, p} = C;
    Cs{2, p} = U;
    Cs{3, p} = diag(Snew);
    Cs{4, p} = diag(S);    
end

%%% blending

%% create morph sequence data
L = sum(pattDurations) + sum(pattTransitions);
mus = zeros(nP,L);

% this initialises the mus in a really weird way
% but it mostly just ends up being an array of ones and zeros
% I'm not really sure what the purpose is...
for window = 1:length(pattOrder)
    if window == 1 % no transition
        mus(pattOrder(1),1:pattDurations(1)) = ...
            ones(1,pattDurations(1));
        startT = pattDurations(1) + 1;
    else
        % start with transition
        mus(pattOrder(window-1),...
            startT:startT+pattTransitions(window-1)-1) = ...
            (pattTransitions(window-1):-1:1) / pattTransitions(window-1);
        mus(pattOrder(window),...
            startT:startT+pattTransitions(window-1)-1) = ...
            (1:pattTransitions(window-1)) / pattTransitions(window-1);
        startT = startT + pattTransitions(window-1);
        mus(pattOrder(window),...
            startT:startT+pattDurations(window)-1) = ...
            ones(1,pattDurations(window));
        startT = startT + pattDurations(window);
    end
end
mus  = smoothmus( mus );

p_CTestPLMorph = zeros(pattDim, L);
x = startXs(:,pattOrder(1));
for n = 1:L
    xOld = x;
    x = (1-LR)*xOld + LR * tanh(W *  x + Wbias);
    % find which mu indices are not 0
    % okay, but why?
    thismu = mus(:,n);
    allInds = (1:nP)';
    muNot0Inds = allInds(thismu ~= 0);
    
    % I think this is used because of different blendy points?
    % and this is operating the blendy points as seperate from the rest
    % of the simulation
    % but I only have a single blendy point
    % also I really don't think I want to use their smoothing paradigm
    if length(muNot0Inds) == 1
        thisC = Cs{1,muNot0Inds};
    else
        thisC = thismu(muNot0Inds(1)) * Cs{1,muNot0Inds(1)} + ...
            thismu(muNot0Inds(2)) * Cs{1,muNot0Inds(2)};
    end
    
    x = thisC * x;
    p_CTestPLMorph(:,n) = Wout * [x; 1];
end

%%% plotting

% the overview on all patterns
figure(1); clf;
set(gcf, 'WindowStyle','normal');
set(gcf,'Position', [800 100 800 800]);

subplotCounter = 0;

for p = 1:nP
    subplotCounter = subplotCounter + 1;
    subplot(nP, 1, subplotCounter);

    hold on;
    % blending results
    plot(1:CtestLength, p_CTestPLSingle(1,:,p)', 'r');
    % after initial training
    %plot(1:size(outsTrain{p},1), outsTrain{p}, 'g');
    hold off;

    % set(gca, 'YLim', [-1 1]);
end