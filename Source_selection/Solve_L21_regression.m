function [M] = Solve_L21_regression(L,W,V,GridSize,IsFine)
  % This function iteratively solve the problem
  %           \argmin_{C_L21} (1/2)||V-L*C_L21*W||_F^2 + alpha*||C_L21_{i}^{T}||_(2,1)   ....(1)
  % where ||x||_(2,1/2) is a composite quasi-norm function.
  % This problem will be solved GridSize times to vary sparsity level by
  % varying parameter alpha in a logspace.
  % input : L, dimension R^{r*m}
  %         W, dimension R^{n*k}
  %         V, dimension R^{r*k}
  %         GridSize: Number of varying regularization parameters between
  %         (10^-5,1)*alpha_max in logspace
  %         IsFine: 1 -> add additional 25 Fine grid alpha
  %               , 0 -> Do nothing
  % alpha_max is the lowest regularization parameter that make the solution
  % of (1) all zero.
  %  output : M.C_L21, parameter C estimated with (2) at all value of alpha [ regularized solution]
  %           M.C_L21_CLS, constrained least-squares solution with zero pattern in parameter C_L21
  %           M.alpha, regularization grid
  %           M.nz_ind_C_L21, nonzero row indices of C_L21
  %           M.ind_chosen_L21, [aicc, bic] selected indices for C_L21_CLS
  % This function solve (1) with spectral_ADMM_sseeg.m
  % Written by PARINTHORN MANOMAISAOWAPAK
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
    [x,~, ~] = spectral_ADMM_sseeg(G, b, P,alpha(ii),0,2,1, [m*n,n,1], 10,0.05,10);
%     [x,~] = grouplasso_sseeg(G, b, P,alpha(ii), n,10);
    C_L21(:,:,ii) = reshape(x,[n m])';
    x0 = x;

else
    [x,~, ~] = spectral_ADMM_sseeg(G, b, P,alpha(ii),0,2,1, [m*n,n,1], 10,0.05,10,x0);
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
