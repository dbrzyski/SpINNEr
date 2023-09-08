% Get estimation error from Simulation Scenario 2
% Input: estimated coefficient matrices from S2_simulated_varyn.m
% Output: One txt file for each method that contains the estimation error
% from each replicate

%addpath('/spinnerCode')
%saddpath('SimulationFunctions/')
%addpath('/glmnet_matlab')


%% Setting
p   = 60;
ns  = [100 150 200 250 300];
nns = length(ns);
ss  = 1;
B1type = 'allS';


%% Compute the corresponding estimation errors.
estErrS = zeros(100, nns); % SPINNER
estErrE = zeros(100, nns); % Elastic Net
estErrN = zeros(100, nns); % Nuclear
estErrL = zeros(100, nns); % Lasso
estErrR = zeros(100, nns); % Ridge
estErrC = zeros(100, nns); % CPM

for i = 1:nns
    if strcmp(B1type,'ones') 
        B1  = ones(8);
    else
        B1  = ss*ones(8);
    end
    B2     = -ss*ones(8);
    B3     = ss*ones(8);
    B      = blkdiag(B1, B2, B3, zeros(p-24, p-24));
    Bnorm  = norm(B - diag(diag(B)), 'fro');
    estErrS(:,i) = estErr_n('S', i, p, B, Bnorm);
    estErrE(:,i) = estErr_n('E', i, p, B, Bnorm);
    estErrN(:,i) = estErr_n('N', i, p, B, Bnorm);
    estErrL(:,i) = estErr_n('L', i, p, B, Bnorm);
    estErrR(:,i) = estErr_n('R', i, p, B, Bnorm);
    estErrC(:,i) = estErr_n('C', i, p, B, Bnorm);
end

dlmwrite('S2_varyn_estErrS.txt', estErrS, '-append', 'delimiter', '\t');
dlmwrite('S2_varyn_estErrE.txt', estErrE, '-append', 'delimiter', '\t');
dlmwrite('S2_varyn_estErrN.txt', estErrN, '-append', 'delimiter', '\t');
dlmwrite('S2_varyn_estErrL.txt', estErrL, '-append', 'delimiter', '\t');
dlmwrite('S2_varyn_estErrR.txt', estErrR, '-append', 'delimiter', '\t');
dlmwrite('S2_varyn_estErrC.txt', estErrC, '-append', 'delimiter', '\t');