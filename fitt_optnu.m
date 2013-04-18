function [f df] = fitt_optnu(nu, delta, p)
% function to solve for ML nu estimate
nu2 = nu/2;
pnu2 = (p+nu)/2;
w = (p + nu) ./ (delta + nu);

f = -psi(nu2) + log(nu2) + (sum(log(w)-w)/length(delta)) + 1 ...
     + psi(pnu2) - log(pnu2);

% jacobian

if nargout > 1
    dp = delta - p;
    dnu = delta + nu;
    sumterm = (1/(length(delta)*(p+nu))) * sum( ((delta-p)./(delta+nu)).^2 );
    df = -0.5*psi(1,nu2) + (1/nu) + sumterm + 0.5*psi(1,pnu2) - 1/(p+nu);
end
