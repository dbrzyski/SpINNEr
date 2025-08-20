function [B_est, beta_est, count] = ProxyLogistic(SVDAx, X, y, Dk, Wk, delta, tol, max_iter, B0, beta0)

% IRLS solver for matrix-on-scalar logistic regression with Frobenius regularization% Inputs:
%   US      - n x n the product of U and S taken from SVD decompositon of A
%   V1      - n x p^2 right eigenvectors of A (only spinning the kernel space of A)
%   X       - n x d matrix of demographic covariates
%   y       - n x 1 binary response vector
%   delta   - regularization parameter
%   tol     - convergence tolerance
%   max_iter - maximum number of IRLS iterations 
%   v1       = V1'*vecY , for Y being a given matrix
%   v2       = (eye(p^2) - V1*V1')*vecY , for Y being a given matrix
%   B0      - (optional) if provided, then B0 is used for initialization of B
%   beta0   - (optional) if provided, then beta0 is used for initialization of beta

% Outputs:
%   B_est   - estimated p x p matrix B
%   beta_est - estimated d x 1 vector beta
%   count   - number of IRLS loops
% -------------------------------------------------------------------------

% create Y 
Y = Dk + Wk/delta;

if SVDAx.UseSymmetricityAndZeroDiag
    delta = 2*delta;
end

% denote SVD matrices
V1 = SVDAx.V;
US = SVDAx.US;
p = SVDAx.pOriginal;
n = size(X,1);
% pp = size(V1,1); % pp = p^2 if UseSymmetricityAndZeroDiag == true, pp = p*(p-1)/2 if UseSymmetricityAndZeroDiag == false

if nargin > 8
    % vecB0 = reshape(B0, [], 1);
    vecB = reshape(B0, [], 1);
    if SVDAx.UseSymmetricityAndZeroDiag
        vecB = vecB(SVDAx.idxs);
    end
    AvecB = US*(V1'*vecB);
else
    vecB = zeros(size(V1,1), 1);
    AvecB = zeros(size(X,1),1);
end

if nargin > 9
    beta = beta0;
else
    beta = zeros(size(X,2), 1);
end

% create two auxilliary vectors
vecY = reshape(Y, [], 1);
diagY = diag(Y);
if SVDAx.UseSymmetricityAndZeroDiag
    vecY = vecY(SVDAx.idxs);
end
v1 = SVDAx.V'*vecY;
V1t_vecY = V1'*vecY;
V1V1t_vecY = V1*V1t_vecY;
v2 = vecY - V1V1t_vecY;
% v2 = (eye(pp) - V1*V1')*vecY;

count = 0;
for iter = 1:max_iter
    eta = AvecB + X*beta;
    mu = 1 ./ (1 + exp(-eta));
    w = mu .* (1 - mu);
    z = eta - (mu - y) ./ (w + eps); % adjusted response
    
    % update z_tilde
    WX = w.*X; % w.*X = W*X
    W2 = diag(w) - WX*((X'*WX)\WX');
    Q1 = US'*W2*US + delta*eye(n);
    z_tilde = Q1\(US'*W2*z + delta*v1);

    % Generates new estimates vectB and beta
    vecB_new = V1*z_tilde + v2;
    AvecB = US*z_tilde;
    beta_new = (X'*WX)\(WX'*(z-AvecB));

    % Check convergence
    if norm(vecB_new - vecB) + norm(beta_new - beta) < tol
        break;
    else
        count = count + 1;
    end
    vecB = vecB_new;
    beta = beta_new;
end

% Reshape and collect B and beta

% Reconstruct the symmetric matrix
if SVDAx.UseSymmetricityAndZeroDiag
    B_est         = zeros(p^2, 1);
    B_est(SVDAx.idxs,:) = vecB; 
    B_est         = reshape(B_est, [p, p]); % p-by-p matrix
    B_est         = B_est + B_est';
    B_est         = B_est + diag(diagY);
else
    B_est = reshape(vecB, p, p);
end
beta_est = beta;

end