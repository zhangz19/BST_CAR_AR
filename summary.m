function [] = summary()

nmod = 4; additive = 0;
datapath = './'; %mat/';
ynames = {'SED load(tons)','SED consentration(mg/kg)','Flowrate(m3/sec)'};
nams = {'Intercept', 'Total_Area', '%Urb', '%Forest', '%Ag', '%%A-soil', '%B-soil',...
        '%C-soil', 'precip(mm)', ...
    '\delta^2', '\tau^2', '\sigma^2', '\gamma', '\phi'};
lambdas = [0.12581, 0.059438, -0.072563];
for nvar = 1:1  %Your Y could be multi-variate
    load('yourData.mat')  % X, Y, W
    Y = Y(:,nvar); allY = Y;
    X = X(:,[1,2,3,5,10:12,18]);
    X = [ones(size(X,1),1), zscore(X)];
    N = size(W,1); T = length(Y)/N; p = size(X,2);
    M = diag(sum(W,1));
    indv = kron(ones(N,1),(1:T)');
    indu = kron((1:N)', ones(T,1));
    npara0 = 4;
    if additive == 1
        npara0 = 5;
    end
    chs = [1,2,3]; nch = length(chs);
    niter = 1e3; burn = 0; thin = 1;
    nsample = (niter-burn)/thin; tot = nch*nsample;
    npara1 = p+npara0; matParas = nan(nch, npara1, nsample);
    
    for nmodel = 1:nmod
        allEs = zeros(1,2); allUs = zeros(tot,N); allVs = zeros(tot,T);
        for ch = 1:nch
            load(strcat(datapath, 'out_',num2str(nvar),'_',num2str(nmodel),'_',num2str(chs(ch)),'.mat'))
            %load(strcat(datapath, 'out_',num2str(nmodel),'_',num2str(chs(ch)),'.mat'))
            n0 = size(matPara, 1);
            allEs = allEs + sum(Es((burn+1):thin:n0,:),1);
            for j = 1:npara1
                matParas(ch, j, :) = matPara((burn+1):thin:n0,j)';
            end
            allUs((ch-1)*nsample + (1:nsample),:) = Us((burn+1):thin:n0,:);
            allVs((ch-1)*nsample + (1:nsample),:) = Vs((burn+1):thin:n0,:);
        end
        allEs = allEs./(nch*nsample);
        R1 = psrf(matParas); boxplot(R1)
        Paras = nan(nch*nsample, npara1);
        for j = 1:npara1
            Paras(:,j) = reshape(matParas(:,j,:),[nch*nsample 1]);
            subplot(4,4,j); hist(Paras(:,j)); xlabel(nams{j}, 'FontSize',15); %ylabel(ynames{nvar}, 'FontSize',15); 
            axis tight
        end
        orient landscape
        print('-painters', '-dpsc2', '-r600', strcat('./tex/stHist',num2str(nvar),'_',num2str(nmodel),'.ps')) %./
        close
        mat = zeros(npara1, 4);
        for i = 1:npara1
            mat(i,1) = mean(Paras(:,i)); mat(i,2) = std(Paras(:,i));
            [lb, ub] = FindHPDset(Paras(:,i), 0.95, []);
            if isempty(lb) || length(lb)>1
                disp('here')
            end
            mat(i,3) = lb(1); mat(i,4) = ub(1);
        end
        disp(mat)
        meanU = mean(allUs, 1); meanV = mean(allVs, 1);
        Yhat = X*mat(1:p,1) + meanU(indu)' + meanV(indv)';
        
        for j = 1:p
            subplot(4,4,j)
            h = plot(X(:,j), Y, 'k.','MarkerSize',17); set(h,'Color',.2*ones(1,3)) %ylim([0,25])
            hold on; h = plot(X(:,j), Yhat, 'r.', 'MarkerSize',10);
            set(h,'Color','r','LineWidth',3); set(gca,'FontSize',15); axis tight
        end
        orient landscape
        print('-painters', '-dpsc2', '-r600', strcat('./tex/stFit',num2str(nvar),'_',num2str(nmodel),'.ps')) %./
        
        RMSE = sqrt(mean((Y - Yhat).^2));
        RSR = sqrt(sum((Y - Yhat).^2))/sqrt(sum((Y - mean(Y)).^2));
        NSE = 1 - sum((Y - Yhat).^2)/sum((Y - mean(Y)).^2);
        Ybar = kron(mean(reshape(Y, [T,N]),1), ones(1,T));
        NSEs = mean(1-sum(reshape((Y - Yhat).^2, [T,N]),1)./sum(reshape((Y - Ybar').^2, [T,N]),1)>0.5)*100;
        PBIAS = sum(Y - Yhat)*100/sum(Y);
        save(strcat(datapath, 'paras',num2str(nvar),'_',num2str(nmodel),'.mat'), ...
            'mat','allEs','RMSE','RSR','NSE','NSEs','PBIAS')
        close
        allY = [allY, Yhat];
    end
    save(strcat(datapath, 'allYs',num2str(nvar),'.mat'), 'allY')
end
end

function [] = summaryOld()

simu = 0; betaprior = 1; logY = 0;
datapath = './'; %../data/LatMonData/';
p = 1; J = 1;
if simu == 1
    J = 1;
end

nch = 3; niter = 5e2; burn = 0; thin = 1;
nsample = (niter-burn)/thin;

for ng = 1:J
    fprintf('group = %d\n', ng)
    load('Ya.mat')
    Y = Y(:,ng); X = X(:,1:4);
    N = size(W,1); T = length(Y)/N; Y = reshape(Y, [T,N]);
    if betaprior == 2
        T = T-1; % make it 2^4
        Y = Y(1:T,:);
        T = size(Y,1); J = floor(log2(T));
        pT = p*T;
        tmp = eye(T); Wav = tmp;
        for i1 = 1:T
            Wav(:,i1) = wavedec(tmp(:,i1),J,'db1');
        end
    elseif betaprior == 1
        pT = size(X,2);
    end
    if logY == 1
        Y = log(Y);
    end
    
    npara0 = 4;
    npara1 = pT+npara0; matParas = nan(nch, npara1, nsample);
    chs = [1,2,3]; nch = length(chs);
    for ch = 1:nch
        load(strcat(datapath, 'out_',num2str((ng-1)*nch+chs(ch)),'.mat'))
        n0 = size(matPara, 1);
        for j = 1:npara1
            matParas(ch, j, :) = matPara((burn+1):thin:n0,j)';
        end
    end
    % R1 = psrf(matParas); boxplot(R1)
    
    Paras = nan(nch*nsample, npara1);
    for j = 1:npara1
        Paras(:,j) = reshape(matParas(:,j,:),[nch*nsample 1]);
    end
    mat = zeros(npara1, 4);
    for i = 1:npara1
        mat(i,1) = mean(Paras(:,i)); mat(i,2) = std(Paras(:,i));
        [lb, ub] = FindHPDset(Paras(:,i), 0.95, []);
        if isempty(lb)
            disp('here')
        end
        mat(i,3) = lb(1); mat(i,4) = ub(1);
    end
    disp(mat(pT+(1:npara0), :))
    
    subplot(1,1,ng); plot(Y, 'Color', [.4,.4,.4])
    hold on;
    if betaprior == 2
        plot(Wav'*mat(1:pT,1), 'Color', [1,0,0]);
    elseif betaprior == 1
        plot(reshape(X*mat(1:pT,1),[T,N]), 'Color', [1,0,0]);
    end
end
end

function [y] = invboxcox(x,lambda)
if lambda == 0 
    y = exp(x);
else
    y = sqrt((lambda*x+1).^(1/lambda - 1));
end
end
