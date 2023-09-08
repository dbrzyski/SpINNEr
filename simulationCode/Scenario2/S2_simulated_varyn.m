%% Simulation scenario 2.
% SPINNER, ElNet, Nuclear, Lasso and Ridge are compared.
% 1) Entries in A_i ~ N(0,1).
% 2) Standardize entries of A_i across subjects, and treat the
% standardized slices as given.
% 3) y_i = <A_i, B> + epsilon_i, where epsilon_i ~ N(0,sigma^2).
% 4) Standardize y, and set X = [].
% B1 = s*ones(8).
% s = 1.
% n varies in 100, 150, 200, 250, 300.

addpath('/spinnerCode')
addpath('SimulationFunctions/')
addpath('/glmnet_matlab')


%% Setting
p   = 60;
ns  = [100 150 200 250 300];
nns = length(ns);
ss  = 1;
B1type = 'allS';


%% Obtain the estimated coefficient matrix from each method.
for i = 1:nns
    fS1 = strcat('S_Bhat_ns', num2str(i), '.txt');
    fS2 = strcat('S_predErr_ns', num2str(i), '.txt');
    fE1 = strcat('E_Bhat_ns', num2str(i), '.txt');
    fE2 = strcat('E_predErr_ns', num2str(i), '.txt');
    fN1 = strcat('N_Bhat_ns', num2str(i), '.txt');
    fN2 = strcat('N_predErr_ns', num2str(i), '.txt');
    fL1 = strcat('L_Bhat_ns', num2str(i), '.txt');
    fL2 = strcat('L_predErr_ns', num2str(i), '.txt');
    fR1 = strcat('R_Bhat_ns', num2str(i), '.txt');
    fR2 = strcat('R_predErr_ns', num2str(i), '.txt');
    fC1 = strcat('C_Bhat_ns', num2str(i), '.txt');
    out = CVall_sim(ns(i), p, B1type, ss, 0.1, 5, fS1, fS2, fE1, fE2, fN1, fN2, fL1, fL2, fR1, fR2, fC1);
end 



