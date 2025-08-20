function out = logistic_spinnerBoth(y, X, SVDAx, lambdaN, lambdaL, WGTs, solOptions)

% This function solves the problem 
%--------------------------------------------------------------------------
%    argmin_{B, beta} {   0.5*sum_i ( y_i - X*beta - <A_i, B> )^2 + 
%                                 lambda_N*|| B ||_* + 
%                               lambda_L*|| vec(B o W) ||_1    }
%--------------------------------------------------------------------------
% for given pair of positive regularization parameters, lambdaN and lambdaL. 
% In order to do so the specific implementation of ADMM algorithm is used.


%% Objects
p = SVDAx.pOriginal;

%% Solver options
deltaInitial1 = solOptions.deltaInitial1;
deltaInitial2 = solOptions.deltaInitial2;
scaleStep     = solOptions.scaleStep;
ratioStep     = solOptions.ratioStep;
mu            = solOptions.mu;
deltaInc      = solOptions.deltaInc;
deltaDecr     = solOptions.deltaDecr;
ratioInc      = solOptions.ratioInc;
ratioDecr     = solOptions.ratioDecr;
maxIters      = solOptions.maxIters;
epsPri        = solOptions.epsPri;
epsDual       = solOptions.epsDual;

%% Initial primal and dual matrix
Dk  = zeros(p,p);
W1k = zeros(p,p);
W2k = zeros(p,p);
B0 = zeros(p,p); 
Beta0 = zeros(size(X,2), 1);

%% ADMM loop
delta1    = deltaInitial1;
delta2    = deltaInitial2;
counterr = 0;
stop     = 0;
CsB      = 0;
DsB      = 0;
DsDp     = 0;
Dlts1    = 0;
Dlts2    = 0;
max_iter = 100;
% tol      = 1e-6;
tol      = 1e-10;

Bnorms = [];
LogProxyCounts = [];
RknormsR1 = [];
RknormsR2 = [];
SknormsR = [];

while stop == 0
    [Bnew, Betanew, logProxyCount] = ProxyLogistic(SVDAx, X, y, Dk, W1k, delta1, tol, max_iter, B0, Beta0);
    Cnew  = ProxG(Dk, W2k, delta2, lambdaN);
    Dnew  = ProxH(Bnew, Cnew, delta1, delta2, W1k, W2k, lambdaL, WGTs);
    W1k   = W1k + delta1*(Dnew - Bnew);
    W2k   = W2k + delta2*(Dnew - Cnew);
    % 
    rk1              = Cnew - Bnew;
    rk2              = Dnew - Bnew;
    sk               = Dnew - Dk;
    rknorm1          = norm(rk1,'fro');
    Bnorm            = norm(Bnew,'fro');
    rknormR1         = rknorm1/Bnorm;
    rknorm2          = norm(rk2,'fro');
    rknormR2         = rknorm2/Bnorm;
    sknorm           = norm(sk,'fro');
    sknormR          = sknorm/norm(Dk,'fro');
    counterr         = counterr + 1;
    CsB(counterr)    = rknormR1; %#ok<*AGROW>
    DsB(counterr)    = rknormR2;
    DsDp(counterr)   = sknormR;
    Dlts1(counterr)  = delta1;
    Dlts2(counterr)  = delta2;
    Dk               = Dnew;  
    B0               = Dnew;
    Beta0            = Betanew;
    
    % ratios update
    if mod(counterr, 20) == 10
        if rknorm1 > mu*rknorm2
            ratioStep = ratioStep*ratioInc;
        else
            if rknorm2 > mu*rknorm1
                ratioStep = ratioStep/ratioDecr;
            end
        end  
    end
    
    % scale update    
    if mod(counterr, 20) == 0
        if mean([rknorm1, rknorm2]) > mu*sknorm
            scaleStep = scaleStep*deltaInc;
        else
            if sknorm > mu*mean([rknorm1, rknorm2])
                scaleStep = scaleStep/deltaDecr;
            end
        end
    end
    delta1 = scaleStep*deltaInitial1;
    delta2 = scaleStep*ratioStep*deltaInitial2;

    % stopping criteria
    if rknormR1 < epsPri && rknormR2 < epsPri && sknormR < epsDual
        stop = 1;
    end
    if counterr > maxIters
        stop = 1;
    end
    if Bnorm <1e-16
        stop = 1;
        Bnew = zeros(p,p);
        Cnew = zeros(p,p);
        Dnew = zeros(p,p);
    end

    Bnorms = [Bnorms, Bnorm];
    LogProxyCounts = [LogProxyCounts, logProxyCount];
    RknormsR1 = [RknormsR1, rknormR1];
    RknormsR2 = [RknormsR2, rknormR2];
    SknormsR = [SknormsR, sknormR];
end

%% Outputs
out         = struct();
out.count   = counterr;
out.Dlts1   = Dlts1;
out.Dlts2   = Dlts2;
out.Blast   = Bnew;
out.Clast   = Cnew;
out.Dlast   = Dnew;
out.B       = Dnew;
out.beta    = Betanew;
out.Bnorms  = Bnorms;
out.RknormsR1  = RknormsR1;
out.RknormsR2  = RknormsR2;
out.SknormsR  = SknormsR;
out.LogProxyCounts  = LogProxyCounts;

end