function [] = getTable3()
% for st
datapath1 = './flow/'; %mat
additive = 0;  % 0 = multiplicative
for nvar = 3:3
    disp(nvar)
    load('Ya.mat')
    Y = Y(:,nvar);
    X = X(:,[1,2,3,5,10:12,18]);
    X = [ones(size(X,1),1), zscore(X)];
    N = size(W,1); T = length(Y)/N; %p = size(X,2);
    % if additive == 1
    %     npara0 = 5;
    % end
    ynames = {'SED_load(tons)','SED-consentration(mg/kg)','Flowrate(m3/sec)'};
    if additive == 1
        nams = {'\\mbox{Intercept}', '\\mbox{total area}', '\\mbox{\\%%Urb}', '\\mbox{\\%%Forest}', ...
            '\\mbox{\\%%Ag}', '\\mbox{\\%%A-soil}', '\\mbox{\\%%B-soil}',...
            '\\mbox{\\%%C-soil}', '\\mbox{precip(mm)}', ...
            '\\delta^2', '\\tau^2', '\\sigma^2', '\\gamma', '\\phi'};
    else
        nams = {'\\mbox{Intercept}', '\\mbox{total area}', '\\mbox{\\%%Urb}', '\\mbox{\\%%Forest}', ...
            '\\mbox{\\%%Ag}', '\\mbox{\\%%A-soil}', '\\mbox{\\%%B-soil}',...
            '\\mbox{\\%%C-soil}', '\\mbox{precip(mm)}', ...
            '\\delta^2', '\\tau^2', '\\phi', '\\gamma'};
    end
    nmod = 4; mats = cell(1,nmod); tmp = zeros(3+5,nmod);
    for model = 1:nmod
        load(strcat(datapath1,'paras',num2str(nvar),'_',num2str(model),'.mat'))
        mats{model} = mat;
        E1 = allEs(1); E2 = allEs(2); DIC = -4*E1 + 2*E2;
        tmp(:,model) = [-2*E1, -2*E1+2*E2, DIC, RMSE, RSR, NSE, NSEs, PBIAS]';
    end
    p = size(mat,1);
    
    for i = 1:p
        fprintf('$'); fprintf(nams{i}); fprintf('$')
        for model = 1:nmod
            if mats{model}(i,1) ~= 0
                fprintf('&$%.2f', mats{model}(i,1)); fprintf('\\,(%.2f, %.2f)$', mats{model}(i,3:4))
            else
                fprintf('&$%.2f$',0)
            end
        end
        fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    end
    fprintf('$\\overline{D(\\utheta)}$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(1,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    fprintf('$p_{D4}$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(2,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    fprintf('$\\mbox{DIC}_4$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(3,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    fprintf('$\\mbox{RMSE}$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(4,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    fprintf('$\\mbox{RSR}$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(5,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    fprintf('$\\mbox{NSE}$'); fprintf(repmat('&$%.4f$',[1,2]), tmp(6,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    fprintf('$\\mbox{NSEs}$'); fprintf(repmat('&$%.4f$',[1,2]), tmp(7,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    fprintf('$\\mbox{PBIAS}$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(8,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
    % disp('done')
end
end

function [] = getTable3old()
% for st
nams = {'\\mu_1', '\\mu_2', '\\mu_3', '\\mu_4', '\\mu_5', '\\mu_6', '\\mu_7', ...
    '\\beta_1', '\\beta_2', '\\beta_3', '\\beta_4', '\\beta_5', '\\beta_6', '\\beta_7',...
    '\\delta^2', '\\tau^2', '\\phi', '\\gamma'};
datapath1 = './st/temporalonly/';
datapath2 = './st/spatialonly/';
datapath3 = './st/nonrandom/';
load(strcat(datapath1,'paras.mat'))
mat1 = mat;
load(strcat(datapath2,'paras.mat'))
mat2 = mat;
load(strcat(datapath2,'paras.mat'))
mat3 = mat;
p = size(mat,1);

for i = 1:p
    fprintf('$'); fprintf(nams{i}); fprintf('$')
    fprintf('&$%.3f', mat1(i,1)); fprintf('\\,(%.3f, %.3f)$', mat1(i,3:4))
    fprintf('&$%.3f', mat2(i,1)); fprintf('\\,(%.3f, %.3f)$', mat2(i,3:4))
    fprintf('&$%.3f', mat3(i,1)); fprintf('\\,(%.3f, %.3f)$', mat3(i,3:4))
    fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
end
tmp = zeros(3,2);
load(strcat(datapath1,'DIC.mat'))
tmp(:,1) = [-2*E1, -2*E1+2*E2, DIC]';
load(strcat(datapath2,'DIC.mat'))
tmp(:,2) = [-2*E1, -2*E1+2*E2, DIC]';
fprintf('$\\overline{D(\\utheta)}$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(1,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
fprintf('$p_{D4}$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(2,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')
fprintf('$\\mbox{DIC}_4$'); fprintf(repmat('&$%.2f$',[1,2]), tmp(3,:)); fprintf('\\\\'); fprintf('\\hline'); fprintf('\n')

disp('done')
end