function [ refSig, outSig ] = varLeaky( patt, LR )

addpath('../helpers');

%%% Experiment control
randstate = 1; % random seed for creating network weights
randn('state', randstate);
rand('twister', randstate);

%%% Setting system params
Netsize = 600; % network size

washoutLength = 50; % how many initial updates per pattern are discarded
% in loading and Wout computation

%%% Ridge regression regularization. Set to 0 because data were so
%%% variable that no regularization was found necessary
TychonovAlpha = 0; % regularizer for  W training
TychonovAlphaReadout = 0; % regularizer for  Wout training

%%% Conceptor learning and testing
alpha = 10; % apertures
CtestLength = 2000;

% weight definitions
Wstar = generate_internal_weights(Netsize, 10/Netsize);
Win =  2*(rand(Netsize, 1) - 0.5);
% WTF? Is this different from what I have in Python? Is my bias totally
% fucked?
Wbias = 2*(rand(Netsize, 1) - 0.5) * 0.8;

learnLength = length(patt) - washoutLength;
xCollector = zeros(Netsize, learnLength);
xOldCollector = zeros(Netsize, learnLength);
WTargetCollector = zeros(Netsize, learnLength);
pCollector = zeros(1, learnLength);
x = zeros(Netsize, 1);
startXs = zeros(Netsize, 1);

for n = 1:(washoutLength + learnLength)
    u = patt(n)'; % pattern input
    xOld = x;
    Wtarget = Wstar * xOld + Win * u;
    x = (1-LR)*xOld + LR * tanh(Wtarget + Wbias);
    if n == washoutLength
        startXs = x;
    end
    if n > washoutLength
        xCollector(:, n - washoutLength ) = x;
        xOldCollector(:, n - washoutLength ) = xOld;
        WTargetCollector(:, n - washoutLength ) = Wtarget;
        pCollector(:, n - washoutLength) = u;

    end
end

%%% compute pattern readout
Wout = (pinv(xCollector * xCollector' + ...
    TychonovAlphaReadout * eye(Netsize)) ...
    * xCollector * pCollector')';
% training error
outsRecovered = Wout * xCollector;
NRMSE_readout = nrmse(outsRecovered, pCollector);
disp(sprintf('mean NRMSE readout: %g   mean abs Wout: %g', ...
    mean(NRMSE_readout(not(isnan(NRMSE_readout)),1)), ...
    mean(mean(abs(Wout)))));

%%% compute W
W = (pinv(xOldCollector * xOldCollector' + ...
    TychonovAlpha * eye(Netsize)) * xOldCollector * WTargetCollector')';
% training errors per neuron
NRMSE_W = nrmse(W * xOldCollector, WTargetCollector);
disp(sprintf('mean NRMSE W: %g   mean abs W: %g ', ...
    mean(NRMSE_W), mean(mean(abs(W)))));

%%% compute conceptors
Cs = cell(4, 1);
pattR = xCollector * xCollector' / learnLength;
[U S V] = svd(pattR);
Snew = (S * pinv(S + alpha^(-2) * eye(Netsize)));
C = U * Snew * U';
Cs{1, 1} = C;
Cs{2, 1} = U;
Cs{3, 1} = diag(Snew);
Cs{4, 1} = diag(S);

%%% test
p_CTestPLSingle = zeros(1, CtestLength);

C = Cs{1, 1};
x = startXs;

for n = 1:CtestLength
    x_rec = x;
    xOld = C * x_rec;
    x = (1-LR)*xOld + LR * tanh(W*xOld + Wbias);

    p_CTestPLSingle(:,n) = Wout * x; 
end

refSig = patt;
outSig = p_CTestPLSingle;

end

