%% Cross validation for elastic-net: tuning both lambda and alpha
%  Note that the objective function for glmnet is
%  0.5 * sum_i (y_i - <A_i, B>)^2/n + 
%  lambda * ( (1-alpha)/2 * ||vec(B)||_2^2 + alpha * ||vec(B)||_1 )

function out = ElNetCV(y, AA, GroupsIdxs)
% GroupsIdxs are the cross-validation indices from SPINNER_CVMESH
% AA and y are assumed to be standardized.

%% Objects
p       = size(AA, 1);
n       = length(y);
Avecs   = reshape(AA, [p^2, n]);
idxs    = logical(reshape(triu(ones(p,p), 1),[p^2,1]));
AvecsUp = 2*Avecs(idxs,:); % (p^2-p)/2-by-n

%% CV options
nalpha  = 15;
alphas  = (0:(nalpha-1))/(nalpha-1);
nlambda = 15;
kfolds  = 5;

%% Cross-validation: use the same indices for different values of alpha
logliksCV = NaN(nlambda, nalpha); % the actual lambda sequences may not be of length nlambda
coefs     = zeros((p^2-p)/2, nalpha);
lambdas   = NaN(nlambda, nalpha);

for ii = 1:nalpha
    opts  = struct('alpha', alphas(ii), 'nlambda', nlambda);
    cvfit = cvglmnet(AvecsUp', y, [], opts, [], kfolds, GroupsIdxs);
    nlam  = length(cvfit.lambda);
    logliksCV(1:nlam,ii) = flip(cvfit.cvm); % mean cross-validated error: 0.5*RSS/n
    coef                 = cvglmnetCoef(cvfit,'lambda_min');
    coefs(:,ii)          = coef(2:end);
    lambdas(1:nlam,ii)   = flip(cvfit.lambda); % the returned lambda sequence is decreasing
end

%% Optimal lambda and alpha
[MM, II]   = min(logliksCV); 
cindex     = find(MM == min(MM)); 
rindex     = II(cindex); 
bestLambda = lambdas(rindex, cindex);
bestAlpha  = alphas(cindex);

%% Final estimate
B         = zeros(p^2, 1);
B(idxs,:) = coefs(:, cindex); 
B         = reshape(B, [p, p]); 
B         = B + B';

%% Output
out              = struct;
out.logliksCV    = logliksCV;
out.lambdas      = lambdas;
out.alphas       = alphas;
out.bestLambda   = bestLambda;
out.bestAlpha    = bestAlpha;
out.B            = B;

end