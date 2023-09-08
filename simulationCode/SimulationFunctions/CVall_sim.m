%% Estimation of coefficient matrix: simulated slices and y
%  Methods considered: spinner, ElNet, nuclear, lasso, ridge, and CPM

%  1) Entries in A_i ~ N(0,1).
%  2) Standardize entries of A_i across subjects, and treat the
%  standardized slices as given.
%  3) y_i = <A_i, B> + epsilon_i, where epsilon_i ~ N(0,sigma^2).
%  4) Standardize y, and set X = [].

function out = CVall_sim(n, p, B1type, s, sigma, nrepl, fS1, fS2, fE1, fE2, fN1, fN2, fL1, fL2, fR1, fR2, fC1)

% B1type: entries in the first block of B; takes on 'ones' or 'allS'
% s: signal strength
% sigma: noise level
% fS1, fS2, fE1, fE2, fN1, fN2, fL1, fL2, fR1, fR2: names of files that
% store the output; 1 for optimal Bhat, 2 for CV estimate of the
%  mean prediction error
% There is no CV step in CPM and p-value threshold is set to 0.01

%% True B
if strcmp(B1type,'ones') 
    B1  = ones(8);
else
    B1  = s*ones(8);
end

B2      = -s*ones(8);
B3      = s*ones(8);
B       = blkdiag(B1, B2, B3, zeros(p-24, p-24));

idxs    = logical(reshape(triu(ones(p,p), 1),[p^2,1])); 
nup     = (p^2-p)/2;

%% SPINNER, elastic-net, nuclear, lasso, ridge
for r = 1:nrepl
    %% Generate upper triangular part of slices  
    SlicesUp     =  randn(nup, n);
    
    %% Standardize slices
    ZSlicesUp    =  zscore(SlicesUp, 0, 2); 
    SlicesVecs   =  zeros(p^2, n);  
    SlicesVecs(idxs,:) = ZSlicesUp;
    Slices       =  reshape(SlicesVecs, [p, p, n]);
    Slices_stdzd =  Slices + permute(Slices, [2 1 3]);
    
    %% Generate y 
    Slices_stdzd_vecs = reshape(Slices_stdzd, [p^2, n]);
    mu                =  Slices_stdzd_vecs' * B(:);
    y                 =  mu + sigma*randn(n,1);

    %% Standardize y
    y_sd         =  std(y);
    y_stdzd      =  zscore(y); 
    
    %% SPINNER, ElNet, nuclear, lasso, ridge
    outSPINNER  =  spinnerCV(y_stdzd, [], Slices_stdzd);
    % Use the GroupsIdxs from SPINNER for all other methods
    GroupsIdxs  =  outSPINNER.GroupsIdxs;
    outElNet    =  ElNetCV(y_stdzd, Slices_stdzd, GroupsIdxs);
    outNuclear  =  spinnerNuclearCV(y_stdzd, [], Slices_stdzd, GroupsIdxs);
    outLasso    =  spinnerLassoCV(y_stdzd, [], Slices_stdzd, GroupsIdxs);
    outRidge    =  ridgeCV(y_stdzd, Slices_stdzd, GroupsIdxs);

    %% CPM, p-value threshold is set to 0.01
    outCPM      =  CPM(y_stdzd, Slices_stdzd, 0.01);
    
    %% Transform Bhats to the original scale
    out = struct;
    out.S_B = y_sd * outSPINNER.B;
    out.E_B = y_sd * outElNet.B;
    out.N_B = y_sd * outNuclear.B;
    out.L_B = y_sd * outLasso.B;
    out.R_B = y_sd * outRidge.B;
    out.C_B = y_sd * outCPM.B;
    
    %% Transform mean cross-validated error (logliksCV) to the original scale
    out.S_logliksCV = y_sd^2 * outSPINNER.logliksCV;
    out.E_logliksCV = y_sd^2 * outElNet.logliksCV;
    out.N_logliksCV = y_sd^2 * outNuclear.logliksCV;
    out.L_logliksCV = y_sd^2 * outLasso.logliksCV;
    out.R_logliksCV = y_sd^2 * outRidge.logliksCV;
    
    %% Write output to files
    dlmwrite(fS1, out.S_B, '-append', 'delimiter', '\t');
    dlmwrite(fS2, out.S_logliksCV, '-append', 'delimiter', '\t');
    dlmwrite(fE1, out.E_B, '-append', 'delimiter', '\t');
    dlmwrite(fE2, out.E_logliksCV, '-append', 'delimiter', '\t');
    dlmwrite(fN1, out.N_B, '-append', 'delimiter', '\t');
    dlmwrite(fN2, out.N_logliksCV, '-append', 'delimiter', '\t');
    dlmwrite(fL1, out.L_B, '-append', 'delimiter', '\t');
    dlmwrite(fL2, out.L_logliksCV, '-append', 'delimiter', '\t');
    dlmwrite(fR1, out.R_B, '-append', 'delimiter', '\t');
    dlmwrite(fR2, out.R_logliksCV, '-append', 'delimiter', '\t');
    dlmwrite(fC1, out.C_B, '-append', 'delimiter', '\t');

end

end




