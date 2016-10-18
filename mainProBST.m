%% Bayesian Spatio-temporal model with separable CAR-AR covariance structure
% * The code is developed for Hamaamin, Y. A., Nejadhashemi, A. P., Zhang, Z., Giri, S. and Woznicki, S. A. (2016).
% "Bayesian Regression and Neuro Fuzzy Methods Reliability Assessment for Estimating Streamflow".
% * The model is: Y = X*beta + Zu + epsilon, with u: NT by 1 spatio-temporal
% random effects with covariance being the kronecker product of
% A: AR(1) model correlation for time series; and
% D: Conditional Auto-regressive(CAR) model for spatial data
% Both assume Markovian structure and have closed-form, sparse inverse. 
% This function also evaluates DIC4 for mixed-effects model for comparison.
% This function can be used for functional regression using spike-and-slab prior for wavelet 
% * Please contact the authors if there are any questions or implementation issues:
% Zhen Zhang, zhangz19@galton.uchicago.edu.

%% Main function for MCMC
function [] = mainProBST(ID)
global N T NT p pT W M X Y eta2 pis J0 L0 facMu indj betaprior J Q nonrandom nontemporal nonspatial
global alphasig invbetasig alphatau invbetatau gap_gamma gap_phi gammas phis lgamma0 lphi0
B = 16e3; burnin = 15e3;
verbose = 0; simu = 0; simudata = 0;  incint = 1; intonly = 0;
betaprior = 1; tranform = 1; computeDIC = 1;
ch = str2double(num2str(ID));
ch0 = ch;
nChain = 3; nmod = 4;
nvar = ceil(ch/(nChain*nmod));
ch = ch - (nvar-1)*nChain*nmod;
nmodel = ceil(ch/nChain);
ch = ch - (nmodel-1)*nChain;
fprintf('nvar = %d, model = %d, chain = %d:\n', [nvar, nmodel, ch]) %partition jobs
switch nmodel
    case 1 % non-random
        nonspatial = 1; nontemporal = 1; nonrandom = 1;
    case 2 % spatial only
        nonspatial = 0; nontemporal = 1; nonrandom = 0;
    case 3 % temporal only
        nonspatial = 1; nontemporal = 0; nonrandom = 0;
    case 4 % full
        nonspatial = 0; nontemporal = 0; nonrandom = 0;
end

load('Ya.mat')
X = X(:,[1,2,3,5,10:12,18]); Y = Y(:,nvar); N = size(W,1);
T = length(Y)/N; NT = N*T; p = size(X,2);
Y = reshape(Y, [T,N]); Y = reshape(Y',[NT,1]);
X0 = nan(NT,p);
for k = 1:p; X0(:,k) = reshape(reshape(X(:,k), [T,N])', [NT,1]); end
X = X0;
M = sum(W,1); M(M==0) = 1; M = diag(M); p = 1;
if tranform == 1; X = zscore(X); end
if incint == 1; X = [ones(NT,1), X]; end
if intonly == 1; X = ones(NT,1); end
J = 1; Q = NT*J;
J0 = floor(log2(T));
indj = repmat((1:NT)',[J,1]);

switch betaprior
    case 1 % flat prior
        p = size(X,2); pT = p*T;
        invXtX = (X'*X)\eye(p); facMu = invXtX*X'; L0 = chol(invXtX, 'lower');
    case 2 % construct wavelet transformation matrix
        J0 = floor(log2(T)); p = 1; NT = N*T; pT = p*T;
        tmp = eye(T); Wav = tmp; for i1 = 1:T; Wav(:,i1) = wavedec(tmp(:,i1),J0,'db1'); end
        X = Wav';
end

% set priors
meansigma2 = 0.01; varsigma2 = 10^2; alphasig = 2+meansigma2^2/varsigma2; invbetasig = meansigma2*(alphasig-1);
meantau2 = 0.01; vartau2 = 10^2; alphatau = 2+meantau2^2/vartau2; invbetatau = meantau2*(alphatau-1);
invM = inv(M); eigs = eig(sqrt(invM)*W*sqrt(invM));
lgamma = max(1/min(eigs),-1); ugamma = 1/max(eigs);
gap_gamma = 1e-2; gammas = (lgamma+gap_gamma):gap_gamma:(ugamma-gap_gamma); len = length(gammas);
lgamma0 = zeros(1,len);
for i = 1:len; lgamma0(i) = 0.5*T*sum(log(eig(M-gammas(i)*W))); end
gap_phi = 1e-2; phis = (-1+gap_phi):gap_phi:(1-gap_phi);
lphi0 = - 0.5*(T-1)*N*log(1-phis.^2);

% set initial values
x.gamma = 0.7337; %0 
x.phi = 0;
x.tau2 = 0.0583; %1e-10;
x.u = normrnd(0,sqrt(x.tau2),[NT,1]); %zeros(NT,1);
switch betaprior
    case 1
        x.beta = facMu*Y;  x.sigma2 = sum((Y - X*x.beta).^2)./Q;
    case 2
        x.beta = zeros(pT,1); x.sigma2 = 0;
        for i = 1:N
            ywave = wavedec(Y(:,i),J0,'db1');
            x.beta = x.beta + ywave;
            x.sigma2 = x.sigma2 + sum((Y(:,i)-ywave).^2);
        end
        x.beta = x.beta./N;  x.sigma2 = x.sigma2./NT;
end
eta2 = 100+zeros(p,T);
pis = 0.5 + zeros(p,J0+1);

if simu == 1 % for simulation
    x.sigma2 = .49; x.tau2 = 4;  x.gamma = .9; V = M - W.*x.gamma; LV = chol(V,'lower'); LV = LV';
    x.phi = 0.8; Bmat = getBmat(x.phi, T, 'mat'); Bmat = Bmat.mat;
    rng('default'); rng(nvar*10);
    myD = x.tau2.*kron(inv(V), Bmat);
    x.u = reshape(mvnrnd(zeros(NT,1), myD),[T,N]);
    if simudata == 1
        eps = normrnd(0,sqrt(x.sigma2),[T,N]);
        Y = repmat(Wav'*x.beta, [1,N])+x.u+eps; x0 = x;
        save('simudatSFC.mat', 'Y','x0')
        plot(Y, 'Color', [.4,.4,.4])
        hold on; plot(Wav'*x.beta, 'Color', [1,0,0]);
    end
    clear('Bmat', 'myD')
    load('simudatSFC.mat')
end
clear('dat','dat0','Wav')

if nonspatial == 1; x.gamma = 0; M = eye(N); end
if nontemporal == 1; x.phi = 0; end
if nonrandom == 1; x.u = zeros(NT,1); x.tau2 = 0; end

V = M - W.*x.gamma; V1 = V\eye(N); LV = chol(V1,'lower'); LV = LV\eye(N);

% MCMC running
matPara = zeros((B-burnin), p+4); Es = [];
Us = zeros((B-burnin), NT);
if computeDIC == 1; Es = zeros(B-burnin, 2); end

rng('default'); rng(ch0*210);
tic
for b = 1:B
    if verbose == 1;  fprintf('%3.2f %3.2f %3.2f %3.2f %3.2f\n', [x.phi, x.gamma, x.sigma2, x.tau2, mean(x.u)]); end
    [x] = updateBeta(x);
    [x,res] = updateSigma2(x);
    if nonrandom ~= 1
        [x,delta,res] = updateTau2(x,LV,res);
        if nontemporal ~= 1; [x] = updatePhi(x, LV); end
        if nonspatial ~= 1; [x,V,LV] = updateGamma(x); end
        [x] = updateU(x,V,delta,res);
    end
    if(b > burnin)
        matPara((b-burnin),:) = [x.beta', x.sigma2, x.tau2, x.phi, x.gamma];
        Us((b-burnin),:) = x.u';
        if computeDIC == 1; Es(b-burnin,:) = getDIC(x,LV); end
    end
end

CPUtime = toc; CPUtime = CPUtime/60;
fprintf('\n%d iterations are done with elapsed time %.2f minutes.\n', B, CPUtime)
nam = strcat('out_',num2str(nvar),'_',num2str(nmodel),'_',num2str(ch),'.mat');
save(nam,'matPara','Us','Es','CPUtime')
end


%% Update fixed-effects
function [x] = updateBeta(x)
global p J T pT N X Y eta2 pis betaprior L0 facMu indj
switch betaprior
    case 1 % flat prior
        Mu = facMu*(Y-x.u(indj));
        x.beta = sqrt(x.sigma2).*(L0*randn([p,1])) + Mu;
    case 2 % wavelet spike-and-slab prior
        for i = 1:p
            t = 1;
            for j = 0:J
                for k = 1:round(2^(j-1))
                    ind0 = false(1,T); ind0(t) = true; inds = false(1, pT); inds(t + (i-1)*T) = true;
                    if i == 1
                        Mu = sum(sum( repmat(X(:,ind0),[1,N]).*...
                            ( Y - x.u - repmat(X(:,~ind0)*x.beta(~ind0),[1,N]) )./x.sigma2 ));
                        Sigma = 1/( N/x.sigma2 + 1/eta2(i,t) );
                        Mu = Sigma*Mu;
                    end
                    logBF = 0.5*Mu^2/Sigma - 0.5*log(eta2(i,t)) + 0.5*log(Sigma);
                    Odds = exp(logBF)*pis(i,j+1)/(1-pis(i,j+1));
                    alpha = 1;
                    if ~isinf(Odds); alpha = Odds/(1+Odds); end
                    u = rand(1);  x.beta(inds) = 0;
                    if u <= alpha; x.beta(inds) = normrnd(Mu, sqrt(Sigma)); end
                    t = t+1;
                end
            end
        end
end
end

%% Update nugget effects (unexplained variation)
function [x,res] = updateSigma2(x)
global Q alphasig invbetasig Y X indj
res = Y - X*x.beta;
azeros = 0.5*Q + alphasig;
bzeros = (0.5*sum((res-x.u(indj)).^2) + invbetasig)^-1;
x.sigma2 = 1./gamrnd(azeros, bzeros);
end

%% Update sill parameter (variation of random effects)
function [x,delta,res] = updateTau2(x,LV,res)
global NT alphatau invbetatau T N J
U = reshape(x.u,[N,T]); U = U';
uMu = fastKronMulti(LV, x.phi, U);
bzeros = (0.5*uMu + invbetatau)^-1;
azeros = 0.5*NT + alphatau;
x.tau2 = 1./gamrnd(azeros, bzeros);
delta = x.tau2/x.sigma2;
res = sum(reshape(res,[NT,J]),2).*delta;
delta = delta*J;
res = reshape(res, [N,T]); res = res';
end

%% Update temporal correlation for AR(1) time series model
function [x] = updatePhi(x,LV)
global phis gap_phi lphi0 T N
U = reshape(x.u,[N,T]); U = U';
loglike = (-0.5/x.tau2).*fastKronMulti(LV, phis, U) + lphi0;
MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
U0 = rand(1);
cump = cumsum([0 P(1:(end-1))]);
i0 = sum(U0 > cump);
x.phi = phis(1);
if i0 > 1; x.phi = phis(i0-1) + gap_phi/P(i0)*(U0-cump(i0)); end
end

%% Update spatial dependence for spatial CAR model
function [x,V,LV] = updateGamma(x)
global T W lgamma0 gammas gap_gamma M N
U = reshape(x.u,[N,T]); U = U'; fac1 = 0;
for t = 1:(T-1)
    if t == 1;  a = 1;  else  a = 1+x.phi^2;  end
    fac1 = fac1 + a.*(U(t,:)*W*U(t,:)') - 2*x.phi.*(U(t,:)*W*U(t+1,:)');
end
fac1 = fac1 + U(T,:)*W*U(T,:)';
loglike = lgamma0 + (0.5/x.tau2)*fac1/(1-x.phi^2)*gammas;
MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
U0 = rand(1);
cump = cumsum([0 P(1:(end-1))]);
i0 = sum(U0 > cump);
x.gamma = gammas(1);
if i0 > 1; x.gamma = gammas(i0-1) + gap_gamma/P(i0)*(U0-cump(i0)); end
V = M - W.*x.gamma;
V1 = V\eye(N); LV = chol(V1,'lower'); LV = LV\eye(N);
end

%% Update high-dimensional spatio-temporal random effects
function [x] = updateU(x,V,delta,res)
global N T NT
ur = zeros(N,T); tmpmu = zeros(N,T); mu = zeros(N,T);
h = (1-x.phi^2)^-1;   L0 = cell(1,T); TL = cell(1,T);  ur0 = normrnd(0, sqrt(x.tau2), [N,T]);

% forward step
for t = 1:T
    if t == 1
        Star = h*V + delta.*eye(N);
        TL{t} = chol(Star, 'lower'); L0{t} = TL{t}\eye(N);
        ur(:,t) = (ur0(:,t)'*L0{t})';
        TL{t} = -x.phi*h*V*L0{t}'; TinvStar = TL{t}*L0{t};
        tmpmu(:,t) = L0{t}*res(t,:)';
    else
        tmpmu(:,t) = res(t,:)' - TL{t-1}*tmpmu(:,t-1); % TL: t-1
        if t>1 && t<T
            Star = (1+x.phi^2)*h*V + delta.*eye(N) + x.phi*h*TinvStar*V; % TinvStar: t-1
        else
            Star = h*V + delta.*eye(N) + x.phi*h*TinvStar*V;
        end
        TL{t} = chol(Star, 'lower'); L0{t} = TL{t}\eye(N);
        tmpmu(:,t) = L0{t}*tmpmu(:,t);
        ur(:,t) = (ur0(:,t)'*L0{t})';
        TL{t} = -x.phi*h*V*L0{t}'; TinvStar = TL{t}*L0{t}; % update TL, TinvStar to t
    end
end

% backward step
Ls = cell(1,T-1);
for t0 = 1:T
    t = T+1-t0;
    if t == T
        mu(:,t) = L0{t}'*tmpmu(:,t);
    else
        for t2 = (t+1):T
            if t2 == t+1; Ls{t} = L0{t+1}; end
            Ls{t2-1} = -Ls{t2-1}*TL{t}*L0{t};
            ur(:,t) = ur(:,t) + (ur0(:,t2)'*Ls{t2-1})';
        end
        mu(:,t) = L0{t}'*(tmpmu(:,t) - TL{t}'*mu(:,t+1));
    end
end
ur = ur + mu;
x.u = reshape(ur, [NT,1]);
end

%% Compute DIC4 for mixed-effects model
function [Es] = getDIC(x0,LV)
global T W M N indj Q NT nonrandom nontemporal nonspatial p Y X
len = 50; x = x0; Es = zeros(1,2);
U = reshape(x.u,[N,T]); U = U';
res = Y - X*x.beta - x.u(indj);
log1 = -0.5*sum((res).^2)/x.sigma2 - 0.5*Q*log(2*pi*x.sigma2);
if nonrandom == 1
    log2 = 0;
else
    lgamma = 0.5*T*sum(log(eig(M-x.gamma*W)));
    lphi = - 0.5*(T-1)*N*log(1-x.phi^2);
    log2 = (-0.5/x.tau2).*fastKronMulti(LV, x.phi, U) + lphi + lgamma - 0.5*NT*log(2*pi*x.tau2);
end
Es(1) = log1+log2;
betabar = zeros(p,1); sigma2bar = 0; tau2bar = 0; phibar = 0; gammabar = 0;
for i = 1:len
    [x] = updateBeta(x); betabar = betabar + x.beta;
    [x,res] = updateSigma2(x); sigma2bar = sigma2bar + x.sigma2;
    if nonrandom ~= 1
        [x,delta,res] = updateTau2(x,LV,res); tau2bar = tau2bar + x.tau2;
        if nontemporal ~= 1; [x] = updatePhi(x, LV); phibar = phibar + x.phi;  end
        if nonspatial ~= 1; x,V,LV] = updateGamma(x); gammabar = gammabar + x.gamma;  end
    end
end
betabar = betabar/len; sigma2bar = sigma2bar/len; tau2bar = tau2bar/len;
phibar = phibar/len; gammabar = gammabar/len;
res = Y - X*betabar - x.u(indj);
log1 = (-0.5/sigma2bar)*sum((res).^2) - 0.5*Q*log(2*pi*sigma2bar);
if nonrandom == 1
    log2 = 0;
else
    lgamma = 0.5*T*sum(log(eig(M-gammabar*W)));
    lphi = - 0.5*(T-1)*N*log(1-phibar^2);
    V = M - W.*gammabar;
    V1 = V\eye(N); LV1 = chol(V1,'lower'); LV1 = LV1\eye(N);
    log2 = (-0.5/tau2bar)*fastKronMulti(LV1, phibar, U) + lphi + lgamma - 0.5*NT*log(2*pi*tau2bar);
end
Es(2) = log1+log2;
end

%% Util function: fast kronecker product of AR(1) and CAR correlation matrices
function [uMu] = fastKronMulti(LD, phi, ur)
% Evaluate uMu = ur'{kron(invD, invB)}ur
% note ur is T by nr matrix
d = length(phi); % phi is 1 by d vector
[T,nr] = size(ur); uMu = zeros(1,d);
h = 1./sqrt(1-phi.^2);
hphi = phi.*h;
if length(size(LD)) < 3
    for s1 = 1:nr
        delta = zeros(T,d);
        for s2 = 1:s1
            if d > 1
                urs = ur(:,s2);
                urnew = repmat(urs, [1 d]);
                urnew(2:T,:) = kron(h, urs(2:T)) - kron(hphi, urs(1:(T-1)));
                delta = delta + urnew.*LD(s1,s2);
            else
                urs = ur(:,s2); urnew = urs;
                urnew(2:T) = (urs(2:T) - phi*urs(1:(T-1)))*h;
                delta = delta + urnew.*LD(s1,s2);
            end
        end
        uMu = uMu + sum((delta).^2, 1);
    end
else
    for s1 = 1:nr
        delta = zeros(T,size(LD,3));
        for s2 = 1:s1
            urs = ur(:,s2); urnew = urs;
            urnew(2:T) = (urs(2:T) - phi*urs(1:(T-1)))*h;
            delta = delta + kron(reshape(LD(s1,s2,:),[1 size(LD,3)]), urnew);
        end
        uMu = uMu + sum((delta).^2, 1);
    end
end
end

%% Util function: obtain AR(1) correlation information (matrix, inverse, cholesky, determinant)
function [out] = getBmat(phi, T, varargin)
% get information of AR(1)-type correlation matrix
mat_ = any(strcmp(varargin, 'mat'));
inv_ = any(strcmp(varargin, 'inv'));
chol_ = any(strcmp(varargin, 'chol'));
det_ = any(strcmp(varargin, 'det'));
time_ = any(strcmp(varargin, 'time'));

t0 = cputime;

if mat_
    out.mat = zeros(T);
    for t1 = 1:T
        for t2 = 1:T
            out.mat(t1,t2) = phi^(abs(t1-t2));
        end
    end
end

if inv_
    out.inv = zeros(T); h = 1/(1-phi^2); hphi = h*phi; hp2 = h*(1+phi^2);
    out.inv(1,1) = h; out.inv(1,2) = -hphi;
    out.inv(T,T) = h; out.inv(T,T-1) = -hphi;
    for t1 = 2:(T-1)
        out.inv(t1,t1) = hp2; 
        out.inv(t1,t1-1) = -hphi; out.inv(t1,t1+1) = -hphi;
    end
end

if chol_
    out.chol = zeros(T); out.chol(:,1) = phi.^(0:(T-1)); h = sqrt(1-phi^2);
    for t1 = 1:(T-1)
        tmp = phi^(t1-1)*h;
        for t2 = 1:(T-t1)
            out.chol(t2+t1,t2+1) = tmp;
        end
    end
end

if det_
    out.det = (1-phi^2)^(T-1);
end

if time_
    disp('The cputime is:'); disp(cputime - t0)
end
end
