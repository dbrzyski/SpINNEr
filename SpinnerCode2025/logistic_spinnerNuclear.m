function out = logistic_spinnerNuclear(y, X, SVDAx, lambda_N, solOptions)

% this function solves the problem 
%
% argmin_B {  0.5*sum_i (y_i - <A_i, B>)^2 + lambda_N || B ||_*  }
% 
% y and AA are after regressing X out

%% Objects
p = SVDAx.pOriginal;

%% Solver options
deltaInitial  = solOptions.deltaInitial1;
mu            = solOptions.mu;
deltaInc      = solOptions.deltaInc;
deltaDecr     = solOptions.deltaDecr;
maxIters      = solOptions.maxIters;
epsPri        = solOptions.epsPri;
epsDual       = solOptions.epsDual;

%% Initial primal and dual matrix
Ck = zeros(p,p);
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
    [Bnew, Betanew, logProxyCount] = ProxyLogistic(SVDAx, X, y, Ck, Wk, delta, tol, max_iter, B0, Beta0);
    Cnew = ProxG(Bnew, -Wk, delta, lambda_N);
    Wk   = Wk + delta*(Cnew - Bnew);  
    rk               = Cnew - Bnew;
    sk               = Cnew - Ck;
    rknorm           = norm(rk,'fro');
    Bnorm            = norm(Bnew,'fro');
    rknormR          = rknorm/Bnorm;
    sknorm           = norm(sk,'fro');
    sknormR          = sknorm/norm(Ck,'fro'); 
    Ck               = Cnew;
    B0               = Cnew;
    Beta0            = Betanew;
    counterr         = counterr + 1;
    
    % stopping criteria
    if counterr > 10
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
        Cnew = zeros(p,p);
    end
end

%% Outputs
out         = struct();
% out.optVal  = optVal;
out.count   = counterr;
out.delta   = delta;
out.Blast   = Bnew;
out.Clast   = Cnew;
out.B       = Cnew;
out.beta    = Betanew;
out.logProxyCount = logProxyCount;

end