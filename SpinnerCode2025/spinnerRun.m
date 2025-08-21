function out = spinnerRun(y, X, AA, lambdaN, lambdaL, W, Params)

% This function solves the problem 
%--------------------------------------------------------------------------
%    argmin_{B, beta} {   0.5*sum_i ( y_i - X*beta - <A_i, B> )^2 + 
%                                 lambda_N*|| B ||_* + 
%                               lambda_L*|| vec(B o W) ||_1    }
%--------------------------------------------------------------------------
% for given pair of regularization parameters, lambdaN and lambdaL. In 
% order to do so the specific implementation of ADMM algorithm is used.

% The SpINNEr toolbox is free software: you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% The SpINNEr toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
%         Authors:    Damian Brzyski and Xixi Hu
%         Date:       June 27, 2018
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%-------------------------------------------
% REQUIRED INPUTS ARGUMENTS:               -
%-------------------------------------------
% y:                 The n-dimensional vector of responses
%~~~~~~~~~~~~~~~~~~~~~
% X:                 matrix of covariates included in the model but not
%                    penalized. Symbol "[]" should be provided to omit it.
%~~~~~~~~~~~~~~~~~~~~~
% AA:                This should be three-dimensional array with dimensions
%                    p,p and n, where ith slice, i.e. A(:,:,i), is a 
%                    symmetric matrix with zeros on the diagonal.
%                    The alternative form, AA = [A1, A2, A3,...,An],is also
%                    supported
%~~~~~~~~~~~~~~~~~~~~~
% lambdaN:           The first regularization parameter
%~~~~~~~~~~~~~~~~~~~~~
% lambdaL:           The second regularization parameter
%
%-------------------------------------------
% OPTIONAL INPUT ARGUMENTS:                -
%-------------------------------------------
% 'W':              Matrix of weights in the SpINNEr optimization problem.
%                   Symmetric p-by-p matrices with nonnegative entries are
%                   possible. The default W has zeros on the diagonal and
%                   ones on the off-diagonal entries
%~~~~~~~~~~~~~~~~~~~~~

%%
%-------------------------------------------
%              OUTPUTS:                    -
%-------------------------------------------
% Outputs can be obtained from the structure class object, 'out'.
%~~~~~~~~~~~~~~~~~~~~~
% B                 SpINNEr estimate of B
%~~~~~~~~~~~~~~~~~~~~~
% beta              SpINNEr estimate of beta
%~~~~~~~~~~~~~~~~~~~~~
% bestLambdaN       The best lambda_N found by cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% bestLambdaL       The best lambda_L found by cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% LambsNgrid        Values of lambda_N considered in cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% LambsLgrid        Values of lambda_L considered in cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% logliksCV         The matrix of prediction errors estimated by
%                   cross-validation for each considered pair of 
%                   tuning parameters (lambda_Ns are in rows)
%-------------------------------------------

%% Checks
% Checking: symmetricity
if ~isequal(AA, permute(AA, [2 1 3]))
    error('Matrices A_i`s must be symmetric')
end

% Checking: zeros on diagonals
if sum(diag(sum(abs(AA),3))) > 0
    error('Matrices A_i`s must have zeros on diagonals')
end

% Checking: zeros on diagonals
if or(isequal(X, zeros(length(y),1) ), isempty(X) ) %checking if X was provided
    XtXXt = zeros(1, length(y));
    X     = zeros(length(y),1);
else
    XtXXt = (X'*X)\X';
end

% Checking: dimmension check
if size(AA,3) ~= length(y)
    error('The third dimension of 3-way tensor containing Ai`s should be the same as the length of y')
end

% Checking: dimmension check
if size(X,1) ~= length(y)
    error('Number of rows in X and the length of y should coincide')
end

%% Objects
p             = size(AA, 1);
n             = length(y);
if nargin == 5
    if lambdaL > 0
        W          = ones(p,p) - eye(p,p);
    else
        W          = ones(p,p);
    end
end

% get rid of X from the optimization problem
H             = eye(n) - X*XtXXt;
AAmatrix      = reshape(AA, [p^2, n])';
AAtilde       = H*AAmatrix;
AAtilde       = reshape(AAtilde', [p, p, n]); % X regressed out
ytilde        = H*y;                          % X regressed out

%% SVD 
% Convert the [p, p, n] array into a (p^2-p)/2-by-n matrix
Avecs     = reshape(AAtilde, [p^2, n]);
idxs      = logical(reshape(triu(ones(p,p), 1),[p^2,1]));
AvecsUp   = 2*Avecs(idxs,:); % (p^2-p)/2-by-n

% Economy-size SVD 
[U, S, V] = svd(AvecsUp, 'econ'); % U is (p^2-p)/2-by-n, S is n-by-n, V is n-by-n.
Sdiag     = diag(S);
Vt        = V';

% SVD objects
SVDAx       = struct;
SVDAx.U     = U;
SVDAx.Sdiag = Sdiag;
SVDAx.Vt    = Vt;
SVDAx.idxs  = idxs;

%% Cases
solverType = double([lambdaN, lambdaL]>0);
solverType = solverType(1) + 2*solverType(2) + 1;

%% Solver
switch solverType
    case 1
        out = minNormEstim(ytilde, SVDAx);
    case 2
        out = spinnerNuclear(ytilde, SVDAx, lambdaN, Params);
    case 3
        out = spinnerLasso(ytilde, SVDAx, lambdaL, W, Params);
    case 4
        out = spinnerBoth(ytilde, SVDAx, lambdaN, lambdaL, W, Params);
end

%% Estimates
estim      = out.B;
beta       = XtXXt*( y - AAmatrix*estim(:) );    

%% Optimial value
DlastVec   = reshape(estim, [p^2,1]);
DlastVecU  = DlastVec(idxs);
MDlast     = U'*DlastVecU.*Sdiag;
optVal_loglik    = 0.5*norm( ytilde - Vt'*MDlast )^2;
optVal_lamNucl   = lambdaN*sum(svd(estim));
optVal_lamLasso  = lambdaL*sum(W(:).*abs(estim(:)));
optVal           = optVal_loglik + optVal_lamNucl + optVal_lamLasso;


%% Outputs
out.optVal  = optVal;
out.optVal_lamNucl  = optVal_lamNucl;
out.optVal_loglik   = optVal_loglik;
out.optVal_lamLasso = optVal_lamLasso;
out.beta = beta;

end