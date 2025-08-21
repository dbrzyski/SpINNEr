function out = logistic_spinnerLasso(y, X, SVDAx, lambdaL, WGTs, solOptions)

% this function solves the problem 
%
% argmin_B {  0.5*sum_i (y_i - <A_i, B>)^2 + lambda_L || B ||_1  }
% 
% y and AA are after regressing X out

%% Objects
p = SVDAx.pOriginal;

%% Solver options
deltaInitial  = solOptions.deltaInitial2;
mu            = solOptions.mu;
deltaInc      = solOptions.deltaInc;
deltaDecr     = solOptions.deltaDecr;
maxIters      = solOptions.maxIters;
epsPri        = solOptions.epsPri;
epsDual       = solOptions.epsDual;

%% Initial primal and dual matrix
Dk = zeros(p,p);
Wk = zeros(p,p);
B0 = zeros(p,p); 
Beta0 = zeros(size(X,2), 1);

%% ADMM loop
delta    = deltaInitial;
max_iter = 100;
tol      = 1e-6;
counterr = 0;
stop     = 0;

while stop == 0
    % Bnew = ProxFsvd(y, SVDAx, Dk, Wk, delta);
    [Bnew, Betanew, logProxyCount] = ProxyLogistic(SVDAx, X, y, Dk, Wk, delta, tol, max_iter, B0, Beta0);
    Dnew = ProxH_lasso(Bnew, delta, Wk, lambdaL, WGTs);
    Wk   = Wk + delta*(Dnew - Bnew);  
    rk               = Dnew - Bnew;
    sk               = Dnew - Dk;
    rknorm           = norm(rk,'fro');
    Bnorm            = norm(Bnew,'fro');
    rknormR          = rknorm/Bnorm;
    sknorm           = norm(sk,'fro');
    sknormR          = sknorm/norm(Dk,'fro'); 
    Dk               = Dnew;
    B0               = Dnew;
    Beta0            = Betanew;
    counterr         = counterr + 1;
    
    % stopping criteria
    if mod(counterr, 10) == 0
        if rknorm > mu*sknorm
            delta = deltaInc*delta;
        else
            if sknorm > mu*rknorm
                delta = delta/deltaDecr;
            end
        end
    end
    if rknormR < epsPri && sknormR < epsDual
        stop = 1;
    end
    if counterr > maxIters
        stop = 1;
    end
    if Bnorm <1e-16
        stop = 1;
        Bnew = zeros(p,p);
        Dnew = zeros(p,p);
    end
end

%% Outputs
out         = struct;
out.count   = counterr;
out.delta   = delta;
out.Blast   = Bnew;
out.Dlast   = Dnew;
out.B       = Dnew;
out.beta    = Betanew;

end