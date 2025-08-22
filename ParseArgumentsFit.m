function [parsedArgs, usingDefaults] = ParseArgumentsFit(providedArgs)

%% Create parser for collecting parameters
p = inputParser;
p.PartialMatching = 0;
p.CaseSensitive = 1;
p.KeepUnmatched = 1;


% defaultW  = ones(p_size, p_size) - eye(p_size); % Default weights
% defaultW = 'Zeros on diagonal, ones on off-diagonal';
defaultUseParallel   = false;
defaultGridLengthN   = 15;
defaultGridLengthL   = 15;
defaultKfolds        = 5;
defaultDisplayStatus = true;

% ScalarNonNeg = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
% zeroOne = @(x) islogical(x) || ismember(x,[0,1]);
% zeroOneInterval = @(x) (x >= 0) && (x <= 1);
% isinteger = @(x) (floor(x)==x);

%% Optional parameters
% addParameter(p, 'W', defaultW, @(x) and(issymmetric(x), all(all(x>=0))));
addParameter(p, 'LambdaN', [], @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
addParameter(p, 'LambdaL', [], @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
addParameter(p, 'Family', 'Gaussian', @(x) ismember(x,{'Gaussian','Binomial'}));
addParameter(p, 'Method', 'CV', @(x) ismember(x,{'CV','Bayesian'}));

% CV options
addParameter(p, 'UseParallel', defaultUseParallel, @islogical );
addParameter(p, 'gridLengthN', defaultGridLengthN, @(x) and(x>0, x == round(x)));
addParameter(p, 'gridLengthL', defaultGridLengthL, @(x) and(x>0, x == round(x)));
addParameter(p, 'gridParameter', 0.75, @(x) x>0);
addParameter(p, 'kfolds', defaultKfolds, @(x) and(x>0, x == round(x)));
addParameter(p, 'displayStatus', defaultDisplayStatus, @islogical );
addParameter(p, 'initLambda', 1, @(x) x>0); % first considered lambdaL and lambdaN
addParameter(p, 'zeroSearchRatio', 100, @(x) x>0); % decides on the speed of increasing lambdas (for finding the zero estimate)
addParameter(p, 'maxLambAcc', 1e-2, @(x) x>0); % precision for finding regularization parameters which give zero estimate for the first time

% ADMM options
addParameter(p, 'UseSymmetricityAndZeroDiag', true, @islogical); % for now implemented only for logistic Spinner
addParameter(p, 'deltaInitial1', 100, @(x) x>0); % the initial "step length" for the update with nuclear norm (i.e. delta1)
addParameter(p, 'deltaInitial2', 100, @(x) x>0); % the initial "step length" for the update with LASSO norm (i.e. delta2)
addParameter(p, 'scaleStep', 1, @(x) x>0); % the initial scale for updated deltas; the scale is changed in repetitions based on the convergence rates
addParameter(p, 'ratioStep', 1, @(x) x>0);  % the initial ratio between updated deltas; the ratio is changed in repetitions based on the convergence rates
addParameter(p, 'mu', 10, @(x) x>0);  % the maximal acceptable ratio between convergence rates to keep deltas without changes in next iteration
addParameter(p, 'deltaInc', 2, @(x) x>0); % the maximal acceptable ratio between convergence rates to keep deltas without changes in next iteration
addParameter(p, 'deltaDecr', 2, @(x) x>0);  % delta is divided by this parameter when the algorithm decides that it should be decreased 
addParameter(p, 'ratioInc', 2, @(x) x>0); % ratio is multiplied by this parameter when the algorithm decides that it should be increased 
addParameter(p, 'ratioDecr', 2, @(x) x>0);  % ratio is divided by this parameter when the algorithm decides that it should be decreased 
addParameter(p, 'maxIters', 50000, @(x) and(x>0, x == round(x)));  % the maximal number of iterations; this is a stopping criterion if the algorithm does not converge
addParameter(p, 'epsPri', 1e-6, @(x) x>0);
addParameter(p, 'epsDual', 1e-6, @(x) x>0);

%% Parsing parameters
parse(p, providedArgs{:});

% checking wheteher there are some unmatched arguments
unMatched = fields(p.Unmatched);
if ~isempty(unMatched)
    warning(['The following arguments names are invalid and were ignored: ' char(join(unMatched, ', '))])
end
allParams = p.Results;
parsedArgs = allParams;
usingDefaults = p.UsingDefaults;

