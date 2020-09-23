function [M] = Solve_L21_regression(L,W,V,GridSize,IsFine)

AA = L'*V*W'; % For alpha_max computation
alpha_max = max(norms(AA,2,2));
G = kron(L,W');
b = V';
b = b(:);
m = size(L,2);
n = size(W,1);
P = eye(m*n);
Lipschitz = norm(W,2)^2*norm(L,2)^2;
STEP_SIZE = 1/Lipschitz;

if IsFine
%     FineGrid = logspace(-8,-5,26)*alpha_max;
FineGrid = logspace(-5,-3,26)*alpha_max;
    FineGrid(end) = [];
%     alpha = [0,FineGrid,logspace(-5,0,GridSize-1)*alpha_max];
    alpha = [FineGrid,logspace(-3,0,GridSize)*alpha_max];

else
%     alpha = [0,logspace(-5,0,GridSize-1)*alpha_max];
    alpha = [logspace(-5,0,GridSize)*alpha_max];
end

for ii=1:length(alpha)
if ii==1
    [x,~, ~] = spectral_ADMM_sseeg_revise(G, b, P,alpha(ii),0,2,1, [m*n,n,1], 10,0.05,10);
%     [x,~] = grouplasso_sseeg(G, b, P,alpha(ii), n,10);
    C_L21(:,:,ii) = reshape(x,[n m])';
    x0 = x;

else
    [x,~, ~] = spectral_ADMM_sseeg_revise(G, b, P,alpha(ii),0,2,1, [m*n,n,1], 10,0.05,10,x0);
%     [x,~] = grouplasso_sseeg(G, b, P,alpha(ii), [2 1 n],0.1,0.5,C_L21(:,:,ii-1));
    C_L21(:,:,ii) = reshape(x,[n m])';
     x0 = x;
end
end

for ii=1:length(alpha)
    fprintf('Model selection progress: %d / %d \n',ii,length(alpha))
nz_ind_C_L21{ii,1} = find(sum(C_L21(:,:,ii).^2,2));
[C_L21_CLS(:,:,ii),~] = constr_LS_eeg(V,L,W,nz_ind_C_L21{ii,1});

[bic_L21(ii),aicc_L21(ii)] = model_criteria(V,L,W,C_L21_CLS(:,:,ii));
[bic_L21_regLS(ii),aicc_L21_regLS(ii)] = model_criteria(V,L,W,C_L21(:,:,ii));

end
[~,ind_chosen_L21.bic] = min(bic_L21);
[~,ind_chosen_L21.aicc] = min(aicc_L21);

[~,ind_chosen_L21.bic_regLS] = min(bic_L21_regLS);
[~,ind_chosen_L21.aicc_regLS] = min(aicc_L21_regLS);


%[C_Lpq_bic,C_L21_bic,C_Lpq,C_L21,bic_ind_Lpq,bic_ind_L21,alpha]
M.C_L21 = C_L21;
M.C_L21_CLS = C_L21_CLS; % Constrained LS L21
M.alpha = alpha;
M.nz_ind_C_L21=nz_ind_C_L21;
M.ind_chosen_L21 = ind_chosen_L21;
end