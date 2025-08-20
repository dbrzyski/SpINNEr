function out = logistic_minNormEstim(y, X, SVDAx)

%% Objects
p = SVDAx.pOriginal;
Dk  = zeros(p,p);
Wk = zeros(p,p);
max_iter = 100;
tol      = 1e-6;

%% Estimate
[B_est, beta_est, count] = ProxyLogistic(SVDAx, X, y, Dk, Wk, 1e-4, tol, max_iter);

%% Output
out            = struct();
out.B          = B_est;
out.beta       = beta_est;
out.logProxyCount = count;

end
