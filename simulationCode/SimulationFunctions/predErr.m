%% Function that calculates the estimation error
%  ||Bhat - B||_F^2/||B||_F^2, ignoring the diagonal terms

function out = predErr(type, i)

file  = strcat(type, '_predErr_ss', num2str(i), '.txt');
cvms  = dlmread(file);
err   = zeros(100,1);

for j = 1:100
    cvm    = cvms(((j-1)*15+1):j*15,:);
    err(j) = min(cvm(:));
end

out = err;

end