function [B_est, beta_est, count] = proxy_logistic(US, V1, X, y, delta, tol, v1, v2, max_iter, B0, beta0)

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

if nargin > 9
    vecB0 = reshape(B0, [], 1);
    vectB = vecB0;
    AvecB = US*(V1'*vecB0);
else
    vectB = zeros(size(V1, 1), 1);
    AvecB = zeros(size(X,1),1);
end

if nargin > 10
    beta = beta0;
else
    beta = zeros(size(X,2), 1);
end

n = size(X,1);
p = sqrt(size(V1, 1));
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
    vectB_new = V1*z_tilde + v2;
    AvecB = US*z_tilde;
    beta_new = (X'*WX)\(WX'*(z-AvecB));

    % Check convergence
    if norm(vectB_new - vectB) + norm(beta_new - beta) < tol
        break;
    else
        count = count + 1;
    end
    vectB = vectB_new;
    beta = beta_new;
end

% Reshape and collect B and beta
B_est = reshape(vectB, p, p);
beta_est = beta;

end