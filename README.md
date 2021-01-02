# **Sp**arsity **I**nducing **N**uclear-**N**orm **E**stimato**r** (SpINNEr)

## The problem which the software addresses
SpINNEr was designed to deal with the estimation of the coefficients matrix under matrix-variate regression model with scalar response. The method works under the assumption that the true solution is well approximable by sparse and low-rank matrix. The immediate application is finding brain connectivities that are associated with a clinical outcome or phenotype.  The proposed framework regresses a (scalar) clinical outcome on matrix-variate predictors which arise in the form of brain connectivity matrices.

## The method
We define the estimate via solution of unconstrained convex optimization problem. Two various penalty terms are added to the objective function:
* **L1 norm** imposed on the vectorized matrix to introduce entry-wise sparsity
* **Nuclear norm** which yields low-rank solution revealing the hidden structure
The detailed description of the method is available on arXiv: https://arxiv.org/abs/2001.11548

## Files structure


## License
The SpINNEr toolbox is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
The SpINNEr toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
