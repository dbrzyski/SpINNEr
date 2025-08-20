function [parsedArgs, usingDefaults] = ParseArguments(providedArgs)

%% Create parser for collecting parameters
p = inputParser;
p.PartialMatching = 0;
p.CaseSensitive = 1;
p.KeepUnmatched = 1;


% defaultW  = ones(p_size, p_size) - eye(p_size); % Default weights
defaultW = 'Zeros on diagonal, ones on off-diagonal';
% ScalarNonNeg = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
% zeroOne = @(x) islogical(x) || ismember(x,[0,1]);
% zeroOneInterval = @(x) (x >= 0) && (x <= 1);
% isinteger = @(x) (floor(x)==x);

%% Optional parameters
addParameter(p, 'Weights', defaultW, @(x) and(issymmetric(x), all(all(x>=0))));

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

