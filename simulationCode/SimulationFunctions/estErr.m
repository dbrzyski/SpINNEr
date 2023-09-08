%% Function that calculates the estimation error
%  ||Bhat - B||_F^2/||B||_F^2, ignoring the diagonal terms

function out = estErr(type, i, p, B, Bnorm)

file  = strcat(type, '_Bhat_ss', num2str(i), '.txt');
Bhats = dlmread(file);
err   = zeros(100,1);

for j = 1:100
    Bhat      = Bhats(((j-1)*p+1):j*p,:);
    BhatB     = Bhat - B;
    BhatBdiag = BhatB - diag(diag(BhatB));
    err(j)    = (norm(BhatBdiag, 'fro')/Bnorm)^2;
end

out = err;

end