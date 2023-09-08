%% Cross validation for nuclear-norm regression: tuning lambda_N

function out = spinnerNuclearCV(y, X, AA, GroupsIdxs)
% GroupsIdxs are the cross-validation indices from spinnerCV
% AA and y are assumed to be standardized.

%% Objects
n           = length(y);

%% CV options
kfolds          = 5;  
gridLengthN     = 15;   
initLambda      = 1;    % first considered lambdaN
zeroSearchRatio = 100;  % decides on the speed of increasing lambdas (for finding the zero estimate)
maxLambAcc      = 1e-1; % precision for finding regularization parameters which give zero estimate for the first time

%% Finding the maximal lambda N
clambdaN  = initLambda;
ValsLambN = zeros(1,1);
counterr1 = 1;
stopp     = 0;

% finding lambda_N for which matrix of zeros is obtained
while stopp == 0
    out = spinner(y, X, AA, clambdaN, 0);
    if norm(out.B, 'fro') < 1e-16
        stopp = 1;
    end
    ValsLambN(counterr1) = clambdaN;
    clambdaN             = zeroSearchRatio*clambdaN;
    counterr1            = counterr1 + 1;
end

% initial interval for maximal lambda N
if length(ValsLambN) == 1
    lamN1 = 0;
    lamN2 = ValsLambN;
else
    lamN1 = ValsLambN(end-1);
    lamN2 = ValsLambN(end);
end

% finding narrow interval for maximal lambda N
stopp      = 0;
counterr2  = 1;
ValsLambNmax = zeros(1,1);
while stopp == 0
    cLamLmaxNew0  = (lamN1 + lamN2)/2;
    outNew0       = spinner(y, X, AA, cLamLmaxNew0, 0);
    if norm(outNew0.B, 'fro') < 1e-16
        lamN2 = cLamLmaxNew0;
    else
        lamN1 = cLamLmaxNew0;
    end
    ValsLambNmax(counterr2) = lamN2;
    counterr2 = counterr2 + 1;
    if abs(lamN2 - lamN1)/lamN2 < maxLambAcc
        stopp = 1;
    end
end

%% Final lambdaN grid
LambsNgrid  = [0, exp( (1:(gridLengthN-1))*log(lamN2+1)/(gridLengthN-1) )- 1]; % grid is defined exponentially in the interval (0, lamN2)

%% Cross - Validation
logliksCV = zeros(gridLengthN, 1);
parfor ii = 1:gridLengthN
    clambdaN = LambsNgrid(ii);
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
        out_CV          =  spinner(y_trening, X_trening, AA_trening, clambdaN, 0);
        normResCV(gg)   =  0.5*norm(y_test - double(ttt(tensor(out_CV.B), tensor(AA_test), 1:2)) - X_test*out_CV.beta)^2;
    end 
    logliksCV(ii) =  sum(normResCV)/n;
    ii
end

%% Optimal lambdaN
[~, II] = min(logliksCV);
bestLambdaN = LambsNgrid(II);

%% Final estimate
outFinal  =  spinner(y, X, AA, bestLambdaN, 0);

%% Output
out              = struct;
out.ValsLambNmax = ValsLambNmax;
out.LambsNgrid   = LambsNgrid;
out.logliksCV    = logliksCV;
out.B            = outFinal.B;
out.bestLambdaN  = bestLambdaN;

end