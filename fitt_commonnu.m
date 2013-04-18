function [grp_mu grp_S nu] = fitt_commonnu(x,y,Ym)
% FITT_COMMONNU(x,y,Ym)
%
% Fit t-distributions with common d.o.f. to grouped data 
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
grp_S = zeros(Nvar, Nvar, Ym);
grp_chS = zeros(Nvar, Nvar, Ym);
grp_Ntrl = zeros(1,Ym);
for yi=1:Ym
    idx = y==(yi-1);
    grpx = x(idx,:);
    grp_x{yi} = grpx;
    grp_mu(:,yi) = mean(grpx);
    S = cov(grpx);
    grp_S(:,:,yi) = S;
    grp_chS(:,:,yi) = chol(S)';
    grp_Ntrl(yi) = size(grpx,1);
end
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
grp_delta = cell(1,Ym);
grp_cx = cell(1,Ym);
while ~converged && (t < maxiter)
    H_old = H;
    t = t+1;
    
    for yi=1:Ym
        % EM update for mu,S for each group independently
        mu = grp_mu(:,yi);
        S = grp_S(:,:,yi);
        chS = grp_chS(:,:,yi);
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
        % ML estimates of mu, S
        mu = sum(bsxfun(@times, thsx, w)) ./ sum(w);
        % centered with updated mean
        cx = bsxfun(@minus, thsx, mu);
        cxw = bsxfun(@times,cx,sqrt(w));
        S = (cxw'*cxw) ./ grp_Ntrl(yi);
        
        % E step again
        chS = chol(S)';
        M = chS\cx';
        delta = sum(M.*M,1)';
        grp_delta{yi} = delta;
        
        grp_S(:,:,yi) = S;
        grp_mu(:,yi) = mu;
        grp_chS(:,:,yi) = chS;
    end
    
    % CM-2 Step
    delta = cell2mat(grp_delta');
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
