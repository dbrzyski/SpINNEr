% Get estimation error from Simulation Scenario 1
% Input: estimated coefficient matrices from S1_simulated_n150ones.m
% Output: One txt file for each method that contains the estimation error
% from each replicate

%addpath('/spinnerCode')
%saddpath('SimulationFunctions/')
%addpath('/glmnet_matlab')


%% Setting
%n      = 150;
p      = 60;
ss     = 2.^(-3:5); % signal strengths; 0.1 is the noise level
nss    = length(ss);
B1type = 'ones';


%% Compute the corresponding estimation errors.
estErrS = zeros(100, nss); % SPINNER
estErrE = zeros(100, nss); % Elastic Net
estErrN = zeros(100, nss); % Nuclear
estErrL = zeros(100, nss); % Lasso
estErrR = zeros(100, nss); % Ridge
estErrC = zeros(100, nss); % CPM

for i = 1:nss
    if strcmp(B1type,'ones') 
        B1  = ones(8);
    else
        B1  = ss(i)*ones(8);
    end
    B2     = -ss(i)*ones(8);
    B3     = ss(i)*ones(8);
    B      = blkdiag(B1, B2, B3, zeros(p-24, p-24));
    Bnorm  = norm(B - diag(diag(B)), 'fro');
    estErrS(:,i) = estErr('S', i, p, B, Bnorm);
    estErrE(:,i) = estErr('E', i, p, B, Bnorm);
    estErrN(:,i) = estErr('N', i, p, B, Bnorm);
    estErrL(:,i) = estErr('L', i, p, B, Bnorm);
    estErrR(:,i) = estErr('R', i, p, B, Bnorm);
    estErrC(:,i) = estErr('C', i, p, B, Bnorm);
end

dlmwrite('S1_n150ones_estErrS.txt', estErrS, '-append', 'delimiter', '\t');
dlmwrite('S1_n150ones_estErrE.txt', estErrE, '-append', 'delimiter', '\t');
dlmwrite('S1_n150ones_estErrN.txt', estErrN, '-append', 'delimiter', '\t');
dlmwrite('S1_n150ones_estErrL.txt', estErrL, '-append', 'delimiter', '\t');
dlmwrite('S1_n150ones_estErrR.txt', estErrR, '-append', 'delimiter', '\t');
dlmwrite('S1_n150ones_estErrC.txt', estErrC, '-append', 'delimiter', '\t');
