%% Simulation scenario 3.
% SPINNER, ElNet, Nuclear, Lasso and Ridge are compared.
%  1) Real slices.
%  2) Standardize entries of A_i across subjects, and treat the
%  standardized slices as given.
%  3) y_i = <A_i, B> + epsilon_i, where epsilon_i ~ N(0,sigma^2).
%  4) Standardize y, and set X = [].

addpath('/spinnerCode')
addpath('SimulationFunctions/')
addpath('/glmnet_matlab')


%% Setting
% Read in real slices.
load('AAYeo.mat'); % AA [148, 148, 100]

ss    = 2.^(-3:5); % signal strengths; 0.1 is the noise level
nss   = length(ss);
B1type = 'ones';


%% Obtain the estimated coefficient matrix from each method.
for i = 7:9
    fS1 = strcat('S_Bhat_ss', num2str(i), '.txt');
    fS2 = strcat('S_predErr_ss', num2str(i), '.txt');
    fE1 = strcat('E_Bhat_ss', num2str(i), '.txt');
    fE2 = strcat('E_predErr_ss', num2str(i), '.txt');
    fN1 = strcat('N_Bhat_ss', num2str(i), '.txt');
    fN2 = strcat('N_predErr_ss', num2str(i), '.txt');
    fL1 = strcat('L_Bhat_ss', num2str(i), '.txt');
    fL2 = strcat('L_predErr_ss', num2str(i), '.txt');
    fR1 = strcat('R_Bhat_ss', num2str(i), '.txt');
    fR2 = strcat('R_predErr_ss', num2str(i), '.txt');
    fC1 = strcat('C_Bhat_ss', num2str(i), '.txt');
    out = CVall_realAA_3blocks(AA, B1type, ss(i), 0.1, 5, fS1, fS2, fE1, fE2, fN1, fN2, fL1, fL2, fR1, fR2, fC1);
end 


