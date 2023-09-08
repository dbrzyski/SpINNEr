% Get estimation error from Simulation Scenario 2
% Input: estimated coefficient matrices from S2_simulated_varyn.m
% Output: One txt file for each method that contains the estimation error
% from each replicate

%addpath('/spinnerCode')
addpath('SimulationFunctions/')
%addpath('/glmnet_matlab')


%% Setting
% Read in real slices.
%load('AAYeo.mat'); % AA [148, 148, 100]

ss    = 2.^(-3:5); % signal strengths; 0.1 is the noise level
nss   = length(ss);
B1type = 'ones';


%% Compute the corresponding estimation errors.
p       = 148;
estErrS = zeros(100, nss); % SPINNER
estErrE = zeros(100, nss); % Elastic Net
estErrN = zeros(100, nss); % Nuclear
estErrL = zeros(100, nss); % Lasso
estErrR = zeros(100, nss); % Ridge
estErrC = zeros(100, nss); % CPM

for i = 1:nss
    B01     = zeros(56); % 1:56
    if strcmp(B1type,'ones') 
        B1  = ones(6); % 57:62
    else
        B1  = ss(i)*ones(6);
    end
    B02     = zeros(6); % 63:68
    B2      = -ss(i)*ones(5); % 69:73
    B03     = zeros(49); % 74:122
    B3      = ss(i)*ones(8); % 123:130
    B04     = zeros(p-130); % 131:p
    B       = blkdiag(B01, B1, B02, B2, B03, B3, B04);
    Bnorm   = norm(B - diag(diag(B)), 'fro');
    estErrS(:,i) = estErr('S', i, p, B, Bnorm);
    estErrE(:,i) = estErr('E', i, p, B, Bnorm);
    estErrN(:,i) = estErr('N', i, p, B, Bnorm);
    estErrL(:,i) = estErr('L', i, p, B, Bnorm);
    estErrR(:,i) = estErr('R', i, p, B, Bnorm);
    estErrC(:,i) = estErr('C', i, p, B, Bnorm);
end

dlmwrite('S3_realAAYeo_estErrS.txt', estErrS, '-append', 'delimiter', '\t');
dlmwrite('S3_realAAYeo_estErrE.txt', estErrE, '-append', 'delimiter', '\t');
dlmwrite('S3_realAAYeo_estErrN.txt', estErrN, '-append', 'delimiter', '\t');
dlmwrite('S3_realAAYeo_estErrL.txt', estErrL, '-append', 'delimiter', '\t');
dlmwrite('S3_realAAYeo_estErrR.txt', estErrR, '-append', 'delimiter', '\t');
dlmwrite('S3_realAAYeo_estErrC.txt', estErrC, '-append', 'delimiter', '\t');