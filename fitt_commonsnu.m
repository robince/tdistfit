function [grp_mu S nu] = fitt_commonsnu(x,y,Ym)
% FITT_COMMONSNU(x,y,Ym)
%
% Fit t-distributions to grouped data with common covariance and d.o.f.  
% using the ECME algorithm (Lui & Rubin, 1995)
%
% Ym - number of groups
% y - group indicators (integers [0, Ym-1])
% 
% C Liu and D B Rubin, (1995) "ML estimation of the t distribution using EM and
% its extensions, ECM and ECME", Statistica Sinica, 5, pp19-39
% http://www3.stat.sinica.edu.tw/statistica/oldpdf/A5n12.pdf
%
if isvector(x)
    x = x(:);
end
Ntrl = size(x,1);
Nvar = size(x,2);
p = Nvar;

% tolerance for entropy
tol = 1e-8;
% Seperate tolerances for S and nu convergence
% S_tol = 1e-6;
% nu_tol = 1e-5;
maxiter = 400;

% Seperate group data and initial conditions
grp_x = cell(1,Ym);
grp_mu = zeros(Nvar, Ym);
grp_cx = cell(Ym,1);
for yi=1:Ym
    idx = y==(yi-1);
    grpx = x(idx,:);
    grp_x{yi} = grpx;
    mu = mean(grpx);
    grp_mu(:,yi) = mu;
    grp_cx{yi} = bsxfun(@minus, grpx, mu);
end
S = cov(cell2mat(grp_cx));
chS = chol(S)';
nu = 10;

t = 1;
converged = false;
H = zeros(1,Ym);

% history of parameters
% theta = zeros(p+numel(S)+2, maxiter);
% theta(:,1) = [mu S(:)' nu H];

% fsolve options
arg = {
'TolFun', 1e-10
'Jacobian', 'on'
% 'DerivativeCheck', 'on'
'Display', 'off'
% 'Algorithm','levenberg-marquardt'
};
arg = arg';
opt = optimset(arg{:});

p2 = p/2;
grp_w = cell(Ym,1);
while ~converged && (t < maxiter)
    H_old = H;
    t = t+1;
    
    for yi=1:Ym
        % EM update for mu for each group independently
        mu = grp_mu(:,yi);
        
        thsx = grp_x{yi};        
        % E step
        % mahalonobis distance with current params        
        cx = bsxfun(@minus, thsx, mu)';
        M = chS\cx;
        % M is the normalised innovation and M(:,i)'*M(:,i) gives the Mahalanobis
        % distance for each x(:,i).
        delta = sum(M.*M,1)';
        w = (p + nu) ./ (delta + nu);
        
        % CM-1 Step
        % ML estimates of mu
        mu = sum(bsxfun(@times, thsx, w)) ./ sum(w);
        
        % centered with updated mean
        cx = bsxfun(@minus, thsx, mu);
        
        grp_w{yi} = w;
        grp_cx{yi} = cx;
        grp_mu(:,yi) = mu;
    end

    % CM-1 Step continued... (S constant across groups)
    cx = cell2mat(grp_cx);
    cxw = bsxfun(@times,cx,sqrt(cell2mat(grp_w)));
    S = (cxw'*cxw) ./ Ntrl;
    
    % E step again
    chS = chol(S)';
    M = chS\cx';
    delta = sum(M.*M,1)';
    
    % CM-2 Step
    optfun = @(v) fitt_optnu(v, delta, p);
    [nu, ~,flag] = fsolve(optfun, nu, opt);
    if flag < 1
        error('fitt:fsolve did not converge')
    end
    
    % convergence detection
    
    % overall difference in parameters
%     theta(:,t) = [mu S(:)' nu H];
%     converged = sum(abs(theta(:,t)-theta(:,t-1))) < tol;

    % don't care about mean for entropy calculation so just check
    % S and nu have converged into a reasonable range
%     converged = (mean(abs(S(:)-S_old(:))) < S_tol) & (abs(nu - nu_old) < nu_tol);
    
    % use entropy as convergence criteria
    nu2 = nu/2;
    nup2 = (nu+p)/2;
    Hnu = log( ((nu*pi).^p2) * beta(p2, nu2) ) - gammaln(p2) ...
               + nup2*(psi(nup2)-psi(nu2));
    for yi=1:Ym
        H(yi) = sum(log(diag(chS))) + Hnu;
    end
    converged = all(abs(H - H_old) < tol);
end
% theta = theta(:,1:t);

if ~converged
    error('fitt:ECME algorithm did not converge (maxiter exceeded)')
end
