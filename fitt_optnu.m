function f = fitt_optnu(nu, delta, p)
% function to solve for ML nu estimate
nu2 = nu/2;
pnu2 = (p+nu)/2;
w = (p + nu) ./ (delta + nu);

f = -psi(nu2) + log(nu2) + (sum(log(w)-w)/length(delta)) + 1 ...
     + psi(pnu2) - log(pnu2);