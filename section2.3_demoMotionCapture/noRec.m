function [ refSig, outSig ] = noRec( patt, LR )

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
    u = patt(n); % pattern input
    xOld = x;
    Wtarget = Win * u;
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

Wout = (pinv(xCollector * xCollector' + ...
    TychonovAlphaReadout * eye(Netsize)) ...
    * xCollector * pCollector')';

outSig = Wout * xCollector;
refSig = pCollector;

end

