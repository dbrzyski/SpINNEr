%% Cross validation for ridge: tuning lambda
%  Note that the objective function for glmnet is
%  0.5 * sum_i (y_i - <A_i, B>)^2/n + 
%  lambda * ( (1-alpha)/2 * ||vec(B)||_2^2 + alpha * ||vec(B)||_1 )

function out = ridgeCV(y, AA, GroupsIdxs)
% GroupsIdxs are the cross-validation indices from spinnerCV
% AA and y are assumed to be standardized.

%% Objects
p       = size(AA, 1);
n       = length(y);
Avecs   = reshape(AA, [p^2, n]);
idxs    = logical(reshape(triu(ones(p,p), 1),[p^2,1]));
AvecsUp = 2*Avecs(idxs,:); 

%% CV options
nlambda = 15;
kfolds  = 5;

%% Cross - Validation
opts      = struct('alpha', 0, 'nlambda', nlambda);
cvfit     = cvglmnet(AvecsUp', y, [], opts, [], kfolds, GroupsIdxs);
lambdas   = flip(cvfit.lambda); % the returned lambda sequence is decreasing
logliksCV = flip(cvfit.cvm);
coef      = cvglmnetCoef(cvfit,'lambda_min'); 
coef      = coef(2:end);

%% Optimal lambda
bestLambda = cvfit.lambda_min;

%% Final estimate
B         = zeros(p^2, 1);
B(idxs,:) = coef; 
B         = reshape(B, [p, p]); 
B         = B + B';

%% Output
out              = struct;
out.logliksCV    = logliksCV;
out.lambdas      = lambdas;
out.bestLambda   = bestLambda;
out.B            = B;

end