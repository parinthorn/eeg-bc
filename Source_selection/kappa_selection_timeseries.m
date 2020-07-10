function [kappa,ind] = kappa_selection_timeseries(data,q,Timepoints,L,alpha,nstate,Cinit)
% This functions compute kappa score for selection of C's row
%        INPUT: EEG data (#channels x #timepoints)
%               q, Group lasso(q=1), non-convex group penalty(q=0.5)
%               Timepoints = timepoints for each chunk
%               L = Leadfield matrix
%               alpha = array of regularization parameter
%               nstate = #state to be determined by subid
%               Cinit = array of initialization for each regularization
%               with dimension (m x n x length(alpha))
%         OUTPUT: kappa = kappa score of each alpha
%                 ind = maximum kappa score index
nchunk = floor(size(data,2)/Timepoints);
r = size(L,1);
y = reshape(data(:,1:nchunk*Timepoints),[r Timepoints nchunk]);
C_kappa = zeros(size(L,2),nchunk,length(alpha));
for ii=1:nchunk
    e = y(:,:,ii);
    [V,W] = extractVW(e,L,[],2*ceil(nstate/r),nstate);
    AA = L'*V*W';
    Lipschitz = norm(L,2)^2*norm(W,2)^2;
    parfor kk = 1:length(alpha)
        fprintf('KAPPA SELECTION-Chunk: %d/%d, alpha %d/%d\n',ii,nchunk,kk,length(alpha))
        [tmp,~] = nmAPG_ss(L,W,V,alpha(kk),q,0.9/Lipschitz,1,Cinit(:,:,ii));
        C_kappa(:,ii,kk) = norms(tmp,2,2);
    end
end
ALL_COMBI = combnk(1:1:nchunk,2);
% kappa = zeros(length(alpha),1);
kappa = zeros(length(alpha),size(ALL_COMBI,1));
for kk = 1:length(alpha)
    parfor ii = 1:size(ALL_COMBI,1)
        C1 = C_kappa(:,ALL_COMBI(ii,1),kk);
        C2 = C_kappa(:,ALL_COMBI(ii,2),kk);
        kappa(kk,ii) = Cohen_kappa(C1,C2)/size(ALL_COMBI,1);
    end
    kappa = sum(kappa,2);
end
[~,ind] = max(kappa);
end