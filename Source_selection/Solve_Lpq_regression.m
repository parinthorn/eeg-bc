function [C_out] = Solve_Lpq_regression(L,W,V,GridSize,IsFine)
% This function iteratively solve the problem
%           \argmin_{C_Lpq} (1/2)||V-L*C_Lpq*W||_F^2 + alpha*||C_Lpq_{i}^{T}||_(2,1/2) ....(1)
%           \argmin_{C_L21} (1/2)||V-L*C_L21*W||_F^2 + alpha*||C_Lpq_{i}^{T}||_(2,1)   ....(2)
% where ||x||_(2,1/2) is a composite quasi-norm function.
% This problem will be solved GridSize times to vary sparsity level by
% varying parameter alpha in a logspace.
% The (
% input : L, dimension R^{r*m}
%         W, dimension R^{n*k}
%         V, dimension R^{r*k}
%         GridSize: Number of varying regularization parameters between
%         (10^-5,1)*alpha_max in logspace
%         IsFine: 1 -> add additional 25 Fine grid alpha \in (10^-8, 10^-5)
%               , 0 -> Do nothing
% alpha_max is the lowest regularization parameter that make the solution
% of (2) all zero (It will sufficiently make (1) zero too but not necessary).
%  output : C_Lpq_bic, parameter C estimated with (1) selected by BIC
%           C_Lpq, parameter C estimated with (1) at all value of alpha
%           C_L21, parameter C estimated with (2) at all value of alpha
%           bic_ind : index of bic selection
AA = L'*V*W'; % For alpha_max computation
alpha_max = max(norms(AA,2,2));
Lipschitz = norm(L,2)^2*norm(W,2)^2;
if IsFine
    FineGrid = logspace(-8,-5,26)*alpha_max;
    FineGrid(end) = [];
    alpha = [0,FineGrid,logspace(-5,0,GridSize-1)*alpha_max];
else
    alpha = [0,logspace(-5,0,GridSize-1)*alpha_max];
end

for ii=1:length(alpha)
if ii==1
    [C_L21(:,:,ii),~] = nmAPG_ss(L,W,V,alpha(ii),1,0.9/Lipschitz,1);
else
    [C_L21(:,:,ii),~] = nmAPG_ss(L,W,V,alpha(ii),1,0.9/Lipschitz,1,C_L21(:,:,ii-1));
end



end
for ii=1:length(alpha)
    fprintf('Estimation progress: %d / %d \n',ii,length(alpha))
    [C_Lpq(:,:,ii),~] = nmAPG_ss(L,W,V,alpha(ii),0.5,0.9/Lipschitz,1,C_L21(:,:,ii));
end
for ii=1:length(alpha)
    fprintf('Model selection progress: %d / %d \n',ii,length(alpha))
[bic_L21(ii),aicc_L21(ii),C_L21_CLS(:,:,ii),flag(ii)] = model_criteria(V,L,C_L21(:,:,ii),W);
[bic_Lpq(ii),aicc_Lpq(ii),C_Lpq_CLS(:,:,ii),flag(ii)] = model_criteria(V,L,C_Lpq(:,:,ii),W);
nz_ind_C_L21{ii,1} = find(C_L21(:,1,ii));% BUG
nz_ind_C_Lpq{ii,1} = find(C_Lpq(:,1,ii)); % BUG
end
[~,bic_ind_Lpq] = min(bic_Lpq);
% C_Lpq_bic = C_Lpq(:,:,bic_ind_Lpq);

[~,bic_ind_L21] = min(bic_L21);
% C_L21_bic = C_L21(:,:,bic_ind_L21);

%[C_Lpq_bic,C_L21_bic,C_Lpq,C_L21,bic_ind_Lpq,bic_ind_L21,alpha]
C_out.C_Lpq = C_Lpq;
% C_out.C_Lpq_bic = C_Lpq_bic;
% C_out.C_L21_bic = C_L21_bic;
C_out.C_L21 = C_L21;
C_out.bic_ind_Lpq = bic_ind_Lpq;
C_out.bic_ind_L21 = bic_ind_L21;
C_out.alpha = alpha;
C_out.C_L21_CLS = C_L21_CLS; % Constrained LS L21
C_out.C_Lpq_CLS = C_Lpq_CLS;% Constrained LS Lpq
C_out.nz_ind_C_L21=nz_ind_C_L21;
C_out.nz_ind_C_Lpq=nz_ind_C_Lpq;
C_out.flag = flag;
end