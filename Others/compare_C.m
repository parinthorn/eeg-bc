function  M = compare_C(ind,ind_true,n,GridSize)
if GridSize>0
[nexp,nreal] = size(ind);
% TP = zeros(nexp,1);
% FP = zeros(nexp,1);
% FN = zeros(nexp,1);
% TN = zeros(nexp,1);
% GridSize = 75;
    TP = zeros(nexp,GridSize);
    FP = zeros(nexp,GridSize);
    FN = zeros(nexp,GridSize);
    TN = zeros(nexp,GridSize);
for ii=1:nexp
    for kk=1:GridSize
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj}{kk}); % true positive

            wrong_ind = setdiff(ind{ii,jj}{kk},ind_true{ii}); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj}{kk}); % false negative
            TPtmp = length(corr_ind);
            FPtmp = length(wrong_ind);
            FNtmp = length(miss_ind);
            npos = TPtmp+FNtmp;
            TNtmp = n-npos-FPtmp;

            TP(ii,kk) = TP(ii,kk)+TPtmp;
            FP(ii,kk) = FP(ii,kk)+FPtmp;
            FN(ii,kk) = FN(ii,kk)+FNtmp;
            TN(ii,kk) = TN(ii,kk)+TNtmp;
        end

    end
end

TP_all = zeros(nexp,GridSize,nreal);
FP_all = zeros(nexp,GridSize,nreal);
FN_all = zeros(nexp,GridSize,nreal);
TN_all = zeros(nexp,GridSize,nreal);
for ii=1:nexp
    for kk=1:GridSize
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj}{kk}); % true positive
            wrong_ind = setdiff(ind{ii,jj}{kk},ind_true{ii}); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj}{kk}); % false negative
            TPtmp = length(corr_ind);
            FPtmp = length(wrong_ind);
            FNtmp = length(miss_ind);
            npos = TPtmp+FNtmp;
            TNtmp = n-npos-FPtmp;
            TP_all(ii,kk,jj) = TPtmp;
            FP_all(ii,kk,jj) = FPtmp;
            FN_all(ii,kk,jj) = FNtmp;
            TN_all(ii,kk,jj) = TNtmp;
        end

    end
end
M.FPR_all = (FP_all)./(TN_all+FP_all);
M.TPR_all = (TP_all)./(TP_all+FN_all);
M.ACC_all = (TP_all+TN_all)./(TN_all+FP_all+FN_all+TP_all);

M.FPR = (FP)./(TN+FP);
M.TPR = (TP)./(FN+TP);
M.ACC = (TP+TN)./(FN+FP+TN+TP);

FP_total = sum(FP,1);
TP_total = sum(TP,1);
FN_total = sum(FN,1);
TN_total = sum(TN,1);

M.FPR_total = (FP_total)./(TN_total+FP_total);
M.TPR_total = (TP_total)./(FN_total+TP_total);
M.ACC_total = (TP_total+TN_total)./(FN_total+FP_total+TN_total+TP_total);
else
    [nexp,nreal] = size(ind);
% TP = zeros(nexp,1);
% FP = zeros(nexp,1);
% FN = zeros(nexp,1);
% TN = zeros(nexp,1);
% GridSize = 75;
    TP = zeros(nexp,1);
    FP = zeros(nexp,1);
    FN = zeros(nexp,1);
    TN = zeros(nexp,1);
for ii=1:nexp
    for kk=1:1
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj}); % true positive
            wrong_ind = setdiff(ind{ii,jj},ind_true{ii}); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj}); % false negative
            TPtmp = length(corr_ind);
            FPtmp = length(wrong_ind);
            FNtmp = length(miss_ind);
            npos = TPtmp+FNtmp;
            TNtmp = n-npos-FPtmp;

            TP(ii,kk) = TP(ii,kk)+TPtmp;
            FP(ii,kk) = FP(ii,kk)+FPtmp;
            FN(ii,kk) = FN(ii,kk)+FNtmp;
            TN(ii,kk) = TN(ii,kk)+TNtmp;
        end

    end
end

TP_all = zeros(nexp,1,nreal);
FP_all = zeros(nexp,1,nreal);
FN_all = zeros(nexp,1,nreal);
TN_all = zeros(nexp,1,nreal);
for ii=1:nexp
    for kk=1:1
        for jj=1:nreal
            corr_ind = intersect(ind_true{ii},ind{ii,jj}); % true positive
            wrong_ind = setdiff(ind{ii,jj},ind_true{ii}); % false positive
            miss_ind = setdiff(ind_true{ii},ind{ii,jj}); % false negative
            TPtmp = length(corr_ind);
            FPtmp = length(wrong_ind);
            FNtmp = length(miss_ind);
            npos = TPtmp+FNtmp;
            TNtmp = n-npos-FPtmp;
            TP_all(ii,kk,jj) = TPtmp;
            FP_all(ii,kk,jj) = FPtmp;
            FN_all(ii,kk,jj) = FNtmp;
            TN_all(ii,kk,jj) = TNtmp;
        end

    end
end
M.FPR_all = squeeze((FP_all)./(TN_all+FP_all));
M.TPR_all = squeeze((TP_all)./(TP_all+FN_all));
M.ACC_all = squeeze((TP_all+TN_all)./(TN_all+FP_all+FN_all+TP_all));

M.FPR = (FP)./(TN+FP);
M.TPR = (TP)./(FN+TP);
M.ACC = (TP+TN)./(FN+FP+TN+TP);

FP_total = sum(FP,1);
TP_total = sum(TP,1);
FN_total = sum(FN,1);
TN_total = sum(TN,1);

M.FPR_total = (FP_total)./(TN_total+FP_total);
M.TPR_total = (TP_total)./(FN_total+TP_total);
M.ACC_total = (TP_total+TN_total)./(FN_total+FP_total+TN_total+TP_total);
end

return;
end
