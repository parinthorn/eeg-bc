function M = compare_F2(ind,ind_true,n)
% This function computes accuracy measures of 'nreal' realization of GC
% matrix in 'nexp' experiments.
% INPUT : ind{i,j} - cell of jth GC matrix's nonzero index on ith
%                    experiment [size (nexp x nreal)] excluding
%                    diagonal indices.
%      ind_true{i} - true GC matrix nonzero index of ith experiment
%                    [size (nexp,1)] excluding diagonal indices
%                n - dimension of GC matrix (nxn), can be vector size
%                    [nexp,1] if each experiment has different size of GC
%                    matrix
% OUTPUT :       M - struct with fields M.TPR, M.FPR, M.ACC,
[nexp,nreal] = size(ind);
TP = zeros(nexp,nreal);
FP = zeros(nexp,nreal);
FN = zeros(nexp,nreal);
TN = zeros(nexp,nreal);

if length(n)>1
    for ii=1:nexp
        npos = length(ind_true{ii});
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj},'rows'); % true positive
            wrong_ind = setdiff(ind{ii,jj},ind_true{ii},'rows'); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj},'rows'); % false negative
            TP(ii,jj) = length(corr_ind);
            FP(ii,jj) = length(wrong_ind);
            FN(ii,jj) = length(miss_ind);
            TN(ii,jj) = n(ii)^2-n(ii)-npos-FP(ii,jj);
        end
    end
else
    for ii=1:nexp
        npos = length(ind_true{ii});
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj},'rows'); % true positive
            wrong_ind = setdiff(ind{ii,jj},ind_true{ii},'rows'); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj},'rows'); % false negative
            TP(ii,jj) = length(corr_ind);
            FP(ii,jj) = length(wrong_ind);
            FN(ii,jj) = length(miss_ind);
            TN(ii,jj) = n^2-n-npos-FP(ii,jj);
        end
    end
end



M.FPR = (FP)./(TN+FP);
M.TPR = (TP)./(FN+TP);
M.ACC = (TP+TN)./(FN+FP+TN+TP);
M.F1  = 2*TP./(2*TP+FP+FN);

M.TP = TP;
M.TN = TN;
M.FP = FP;
M.FN = FN;
end
