%-------------------------------------------------------------------------- 
% This function finds the CPM estimate for given threshold parameter 
% pv_thresh. The estimate is a result of a following procedure:
% * find p^2 correlations between edges and the response variable
% * locate the edges for which correlation p-values is smaller than 
%   pv_thresh
% * define the final estimate as p x p matrix with zeros corresponding to 
%   indices where p-values > pv_thresh and simple linear estimates for the 
%   idices where p-values < pv_thresh
%  
%  For more details on CPM algorithm see
%  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5526681/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
%         Authors:    Damian Brzyski
%         Date:       September 6, 2021
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%-------------------------------------------
% REQUIRED INPUTS ARGUMENTS:               -
%-------------------------------------------
% y:                 The n-dimensional vector of responses
%
%~~~~~~~~~~~~~~~~~~~~~
% AA:                This should be three-dimensional array with dimensions
%                    p,p and n, where ith slice, i.e. A(:,:,i), is a 
%                    symmetric matrix with zeros on the diagonal.
%                    The alternative form, AA = [A1, A2, A3,...,An],is also
%                    supported
%
%~~~~~~~~~~~~~~~~~~~~~
% pv_thresh:         p-value threshold
%-------------------------------------------

%%
%-------------------------------------------
%              OUTPUTS:                    -
%-------------------------------------------
% Outputs can be obtained from the structure class object, 'out'.
%~~~~~~~~~~~~~~~~~~~~~
% B                 CPM estimate of B
%-------------------------------------------

function out = CPM(y, AA, pv_thresh)

%% Objects
[p, p2, p3] = size(AA);
n = length(y);

%% Tensor of regressor matrices
if p3>1
    AA = reshape(AA, [p^2, n]);
else
    if mod(p2,p)~=0
        error('The assumed form is AA = [A1, A2,...,An], with square matrices Ai, hence number of columns in AA should be the multiplicity of the number of rows')
    else
        AA = reshape(AA, [p, p, p2/p]);
        AA = reshape(AA, [p^2, n]);
    end
end

%% Get the estimate
%AA2 = zscore(AA')';             % standardize data
% AA is assumed to be already standardized in simulations
[estim, pvals] = corr(AA',y);
estim = estim*std(y)/2;  % division by 2 to account for symmetricity 
mask = pvals < pv_thresh;
estim(mask == 0) = 0;
Bhat = reshape(estim, [p, p]);

%% Output
out = struct;
out.B = Bhat;
out.intercept = mean(y);

end


