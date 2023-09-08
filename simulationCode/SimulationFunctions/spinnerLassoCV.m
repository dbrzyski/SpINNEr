%% Cross validation for lasso: tuning lambda_L

function out = spinnerLassoCV(y, X, AA, GroupsIdxs)
% GroupsIdxs are the cross-validation indices from spinnerCV
% AA and y are assumed to be standardized.

%% Objects
n           = length(y);

%% CV options
kfolds          = 5;
gridLengthL     = 15;     
initLambda      = 1;    % first considered lambdaN
zeroSearchRatio = 100;  % decides on the speed of increasing lambdas (for finding the zero estimate)
maxLambAcc      = 1e-1; % precision for finding regularization parameters which give zero estimate for the first time

%% Finding the maximal lambda L
clambdaL  = initLambda;
ValsLambL = zeros(1,1);
counterr1 = 1;
stopp     = 0;

% finding lambda_L for which matrix of zeros is obtained
while stopp == 0
    out = spinner(y, X, AA, 0, clambdaL);
    if norm(out.B, 'fro') < 1e-16
        stopp = 1;
    end
    ValsLambL(counterr1) = clambdaL;
    clambdaL             = zeroSearchRatio*clambdaL;
    counterr1            = counterr1 + 1;
end

% initial interval for maximal lambda L
if length(ValsLambL) == 1
    lamL1 = 0;
    lamL2 = ValsLambL;
else
    lamL1 = ValsLambL(end-1);
    lamL2 = ValsLambL(end);
end

% finding narrow interval for maximal lambda L
stopp      = 0;
counterr2  = 1;
ValsLambLmax = zeros(1,1);
while stopp == 0
    cLamLmaxNew0  = (lamL1 + lamL2)/2;
    outNew0       = spinner(y, X, AA, 0, cLamLmaxNew0);
    if norm(outNew0.B, 'fro') < 1e-16
        lamL2 = cLamLmaxNew0;
    else
        lamL1 = cLamLmaxNew0;
    end
    ValsLambLmax(counterr2) = lamL2;
    counterr2 = counterr2 + 1;
    if abs(lamL2 - lamL1)/lamL2 < maxLambAcc
        stopp = 1;
    end
end

%% Final lambdaL grid
LambsLgrid  = [0, exp( (1:(gridLengthL-1))*log(lamL2+1)/(gridLengthL-1) )- 1]; % grid is defined exponentially in the interval (0, lamL2)

%% Cross - Validation
logliksCV = zeros(gridLengthL, 1);
parfor ii = 1:gridLengthL
    clambdaL = LambsLgrid(ii);
    normResCV = zeros(kfolds,1);
    for gg = 1:kfolds
        testIndices     =  find(GroupsIdxs == gg);
        treningIndices  =  setdiff(1:n , testIndices)';
        AA_trening      =  AA(:,:, treningIndices);
        AA_test         =  AA(:,:, testIndices);
        if ~isempty(X)
            X_trening       =  X(treningIndices,:);
            X_test          =  X(testIndices,:);
        else
            X_trening       =  [];
            X_test          =  0;
        end
        y_trening       =  y(treningIndices);
        y_test          =  y(testIndices);
        out_CV          =  spinner(y_trening, X_trening, AA_trening, 0, clambdaL);
        normResCV(gg)   =  0.5*norm(y_test - double(ttt(tensor(out_CV.B), tensor(AA_test), 1:2)) - X_test*out_CV.beta)^2;
    end 
    logliksCV(ii) =  sum(normResCV)/n;
    ii
end

%% Optimal lambdaN
[~, II] = min(logliksCV);
bestLambdaL = LambsLgrid(II);

%% Final estimate
outFinal  =  spinner(y, X, AA, 0, bestLambdaL);

%% Output
out              = struct;
out.ValsLambLmax = ValsLambLmax;
out.LambsLgrid   = LambsLgrid;
out.logliksCV    = logliksCV;
out.B            = outFinal.B;
out.bestLambdaL  = bestLambdaL;

end