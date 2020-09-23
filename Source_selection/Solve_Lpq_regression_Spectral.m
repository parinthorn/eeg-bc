function [M] = Solve_Lpq_regression_Spectral(L,W,V,GridSize,IsFine)

AA = L'*V*W'; % For alpha_max computation
alpha_max = max(norms(AA,2,2));
G = kron(L,W');
b = V';
b = b(:);
m = size(L,2);
n = size(W,1);
P = eye(m*n);
STEP_SIZE = 1/norm(G,2)^2;

if IsFine
%     FineGrid = logspace(-8,-5,26)*alpha_max;
%     FineGrid = logspace(-3,-2,26)*alpha_max;
FineGrid = logspace(-5,-4,26)*alpha_max;
    FineGrid(end) = [];
%     alpha = [0,FineGrid,logspace(-5,0,GridSize-1)*alpha_max];
%     alpha = [0,FineGrid,logspace(-2,-1,GridSize-1)*alpha_max];
    alpha = [0,FineGrid,logspace(-4,-1,GridSize-1)*alpha_max];

else
    alpha = [0,logspace(-5,0,GridSize-1)*alpha_max];
end

for ii=1:length(alpha)
    d = 0.8;
    eta = 0.7;
    rho = 0.5;
if ii==1
    [x,~, ~] = spectral_ADMM_sseeg(G, b, P,alpha(ii),0,2,1, [m*n,n,1], 10,0.5,2);
    C_L21(:,:,ii) = reshape(x,[n m])';
    x0 = x;
    [C_Lpq(:,:,ii),~] = nmAPG_ss(L,W,V,alpha(ii),0.5,STEP_SIZE,1,C_L21(:,:,ii));

else
    [x,~, ~] = spectral_ADMM_sseeg(G, b, P,alpha(ii),0,2,1, [m*n,n,1], 10,0.5,2,x0);
    C_L21(:,:,ii) = reshape(x,[n m])';
     x0 = x;
    [C_Lpq(:,:,ii),~] = nmAPG_ss(L,W,V,alpha(ii),0.5,STEP_SIZE,1,C_L21(:,:,ii));
end
end


for ii=1:length(alpha)
    fprintf('Model selection progress: %d / %d \n',ii,length(alpha))
nz_ind_C_L21{ii,1} = find(sum(C_L21(:,:,ii).^2,2));
nz_ind_C_Lpq{ii,1} = find(sum(C_Lpq(:,:,ii).^2,2));
[C_L21_CLS(:,:,ii),~] = constr_LS_eeg(V,L,W,nz_ind_C_L21{ii,1});
[C_Lpq_CLS(:,:,ii),~] = constr_LS_eeg(V,L,W,nz_ind_C_Lpq{ii,1});

[bic_L21(ii),aicc_L21(ii)] = model_criteria(V,L,W,C_L21_CLS(:,:,ii));
[bic_Lpq(ii),aicc_Lpq(ii)] = model_criteria(V,L,W,C_Lpq_CLS(:,:,ii));
[bic_L21_regLS(ii),aicc_L21_regLS(ii)] = model_criteria(V,L,W,C_L21(:,:,ii));
[bic_Lpq_regLS(ii),aicc_Lpq_regLS(ii)] = model_criteria(V,L,W,C_Lpq(:,:,ii));

end

[~,ind_chosen_Lpq.bic] = min(bic_Lpq);
[~,ind_chosen_L21.bic] = min(bic_L21);
[~,ind_chosen_Lpq.aicc] = min(aicc_Lpq);
[~,ind_chosen_L21.aicc] = min(aicc_L21);

[~,ind_chosen_Lpq.bic_regLS] = min(bic_Lpq_regLS);
[~,ind_chosen_L21.bic_regLS] = min(bic_L21_regLS);
[~,ind_chosen_Lpq.aicc_regLS] = min(aicc_Lpq_regLS);
[~,ind_chosen_L21.aicc_regLS] = min(aicc_L21_regLS);

M.C_L21 = C_L21;
M.C_Lpq = C_Lpq;
M.C_L21_CLS = C_L21_CLS; % Constrained LS L21
M.C_Lpq_CLS = C_Lpq_CLS;% Constrained LS Lpq
M.alpha = alpha;
M.nz_ind_C_L21=nz_ind_C_L21;
M.nz_ind_C_Lpq=nz_ind_C_Lpq;
M.ind_chosen_Lpq = ind_chosen_Lpq;
M.ind_chosen_L21 = ind_chosen_L21;
end
