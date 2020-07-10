function M = compare_F(ind,ind_true,n,singleROC)
% This function computes accuracy measures of 'nreal' realization of GC
% matrix in 'nexp' experiments.
% INPUT : ind{i,j} - cell of jth GC matrix's nonzero index on ith
%                      experiment [size (nexp x nreal)].
%      ind_true{i} - true GC matrix nonzero index of ith experiment
%                    [size (nexp,1)]
%                n - dimension of GC matrix (nxn), can be vector size
%                   [nexp,1] if each experiment has different size of GC
%                   matrix
% OUTPUT :       M - struct with fields M.TPR, M.FPR, M.ACC,
[nexp,nreal] = size(ind);
TP = zeros(nexp,1);
FP = zeros(nexp,1);
FN = zeros(nexp,1);
TN = zeros(nexp,1);

if nargin<4
    singleROC = 0;
end

if singleROC==1
    for ii=1:nexp
        for jj=1:nreal
            corr_ind = intersect(ind_true,ind{ii,jj},'rows'); % true positive
            wrong_ind = setdiff(ind{ii,jj},ind_true,'rows'); % false positive
            miss_ind = setdiff(ind_true,ind{ii,jj},'rows'); % false negative
            TPtmp = length(corr_ind);
            FPtmp = length(wrong_ind);
            FNtmp = length(miss_ind);
            npos = TPtmp+FNtmp;
            TNtmp = n^2-npos-FPtmp;
            
            TP(ii) = TP(ii)+TPtmp;
            FP(ii) = FP(ii)+FPtmp;
            FN(ii) = FN(ii)+FNtmp;
            TN(ii) = TN(ii)+TNtmp;
        end
    end
    M.FPR = (FP)./(TN+FP);
    M.TPR = (TP)./(FN+TP);
    M.ACC = (TP+TN)./(FN+FP+TN+TP);
    return;
end

if length(n)==1
    for ii=1:nexp
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj},'rows'); % true positive
            wrong_ind = setdiff(ind{ii,jj},ind_true{ii},'rows'); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj},'rows'); % false negative
            TPtmp = length(corr_ind);
            FPtmp = length(wrong_ind);
            FNtmp = length(miss_ind);
            npos = TPtmp+FNtmp;
            TNtmp = n^2-npos-FPtmp;
            
            TP(ii) = TP(ii)+TPtmp;
            FP(ii) = FP(ii)+FPtmp;
            FN(ii) = FN(ii)+FNtmp;
            TN(ii) = TN(ii)+TNtmp;
        end
    end
elseif (length(n)==nexp)
    for ii=1:nexp
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj},'rows'); % true positive
            wrong_ind = setdiff(ind{ii,jj},ind_true{ii},'rows'); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj},'rows'); % false negative
            TPtmp = length(corr_ind);
            FPtmp = length(wrong_ind);
            FNtmp = length(miss_ind);
            npos = TPtmp+FNtmp;
            TNtmp = n(ii)^2-npos-FPtmp;
            
            
            TP(ii) = TP(ii)+TPtmp;
            FP(ii) = FP(ii)+FPtmp;
            FN(ii) = FN(ii)+FNtmp;
            TN(ii) = TN(ii)+TNtmp;
            
        end
    end
else
    error('The dimension n is not valid.')
    
end
M.FPR = (FP)./(TN+FP);
M.TPR = (TP)./(FN+TP);
M.ACC = (TP+TN)./(FN+FP+TN+TP);
return;
end