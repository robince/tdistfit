% FITT(x)
%
% Fit a t-distribution using the ECME algorithm (Lui & Rubin)
%
if isvector(x)
    x = x(:);
end
Ntrl = size(x,1);
Nvar = size(x,2);
p = Nvar;

S_tol = 1e-6;
nu_tol = 1e-5;
tol = 1e-5;
maxiter = 100;
converged = false;

% Initial conditions
mu = mean(x);
S = cov(x);
nu = 10;

t = 1;
H = 0;
% avoid computing full likelihood
% check for convergence of parameters
theta = zeros(p+numel(S)+2, maxiter);
theta(:,1) = [mu S(:)' nu H];
opt = optimset('TolFun',1e-10 ...
               ,'Display','off' ...
               );%,'Algorithm','levenberg-marquardt');

p2 = p/2;
while ~converged && (t < maxiter)
    S_old = S;
    nu_old = nu;
    H_old = H;
    
    % E step
    % mahalonobis distance with current params
    chS = chol(S)';
    cx = bsxfun(@minus, x, mu)';
    M = chS\cx;
    % M is the normalised innovation and M(:,i)'*M(:,i) gives the Mahalanobis
    % distance for each x(:,i).
    delta = sum(M.*M,1)';
    w = (p + nu) ./ (delta + nu);
    
    % CM-1 Step
    % ML estimates of mu, S
    mu = sum(bsxfun(@times, x, w)) ./ sum(w);
    % centered with updated mean
    cx = bsxfun(@minus, x, mu);
    cxw = bsxfun(@times,cx,sqrt(w));
    S = (cxw'*cxw) ./ Ntrl;
    chS = chol(S)';
    
    % line search is slow so only do it every other iteration
    if mod(t+1,2)==0
        t;
        % E step again        
        M = chS\(cx');
        delta = sum(M.*M,1)';
        w = (p + nu) ./ (delta + nu);
        
        % CM-2 Step
        optfun = @(v) fitt_optnu(v, delta, p);
%             nu = fminsearch(optfun,nu,opt);
        nu = fsolve(optfun, nu, opt);
    end
    
    %TODO: proper log likelihood for convergence?
    t = t+1;
    
%     converged = sum(abs(theta(:,t)-theta(:,t-1))) < tol;
    % don't care about mean for entropy calculation so just check
    % S and nu have converged into a reasonable range
%     converged = (mean(abs(S(:)-S_old(:))) < S_tol) & (abs(nu - nu_old) < nu_tol);
    % use entropy as convergence criteria
    nu2 = nu/2;
    nup2 = (nu+p)/2;
    H = sum(log(diag(chS))) ...
               + log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
               + nup2*(psi(nup2)-psi(nu2));
    theta(:,t) = [mu S(:)' nu H];
    converged = abs(H - H_old) < tol;
end
theta = theta(:,1:t);
