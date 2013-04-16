function [c, S, nu] = fitt_approx(x)
% FITT_APPROX(x)
%
% Fit a multivariate t-distribution to data using the approximation method
% from:
% C Aeschlimna, J Park and KA Cak, "A Novel Parameter Estimation Algorithm
% for the Multivariate t-Distribution and its Application to Computer
% Vision" ECCV 2010 
% http://link.springer.com/chapter/10.1007%2F978-3-642-15552-9_43

Ntrl = size(x,1);
Nvar = size(x,2);

c = median(x);
% centered data
cx = bsxfun(@minus, x, median(c));

zi = log(sum(cx.^2,2));
z = sum( zi ) ./ Ntrl;

b = (sum( (zi-z).^2 ) ./ Ntrl) - psi(1, Nvar/2);
nu = (1 + sqrt(1+4*b)) / b;

alpha = exp(z - log(nu) + psi(0, nu/2) - psi(0, Nvar/2));
beta = (2*log2(Nvar))/(nu^2 + log2(Nvar));

S = (cx'*bsxfun(@rdivide,cx, sum(cx.^2,2).^(beta/2)))./Ntrl;
S = (alpha*Nvar/trace(S)) * S;