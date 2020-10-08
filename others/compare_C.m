function  M = compare_C(ind,ind_true,m,GridSize)
% This function performs sparsity comparison of nonzero row 
% indices "ind" and the estimated matrix C with nonzero row indices "ind_true" 
% of ground-truth matrix C.
% This function compare sparsity pattern of "nexp" ground-truth C matrices with 
% "nreal" estimated C matrices. "nreal" estimated C matrices share the same ground-truth C. 
% The sparsity pattern is varied for "GridSize" times.
% If GridSize is set to -1, the single comparison will be made. This is for
% accuracy evaluation of a single C using some criteria such as BIC or AICC
%
% input: ind      - non-zero row linear index of estimated C stored in cell
%                   size (nexp,nreal) where nreal is the number of
%                   instances that share same ground-truth models.
%        ind_true - non-zero row linear index of ground-truth C stored in 
%                   cell size (nexp,1).
%        m        - # row of estimated C
%        GridSize - number of regularization level,
%                   For single regularization level use GridSize = -1
% output: M.FPR_all
% M.FPR_all : cell array of FPR of all realization from all ground-truth with dimension (nexp x GridSize x nreal), 
% M.TPR_all : cell array of TPR of all realization from all ground-truth with dimension (nexp x GridSize x nreal), 
% M.ACC_all : cell array of ACC of all realization from all ground-truth with dimension (nexp x GridSize x nreal), 
% 
% M.FPR = cell array of FPR of each ground-truth in all realization with dimension (nexp x GridSize), 
% M.TPR = cell array of TPR of each ground-truth in all realization with dimension (nexp x GridSize), 
% M.ACC = cell array of ACC of each ground-truth in all realization with dimension (nexp x GridSize), 
% 
% M.FPR_total = cell array of FPR of each realization across all ground-truth with dimension (nreal x GridSize), 
% M.TPR_total = cell array of FPR of each realization across all ground-truth with dimension (nreal x GridSize), 
% M.ACC_total = cell array of FPR of each realization across all ground-truth with dimension (nreal x GridSize), 

%
% Written by PARINTHORN MANOMAISAOWAPAK, 2020
%
% Licenses
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
if GridSize>0
[nexp,nreal] = size(ind);
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
            TNtmp = m-npos-FPtmp;
            
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
            TNtmp = m-npos-FPtmp;
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
            TNtmp = m-npos-FPtmp;
            
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
            TNtmp = m-npos-FPtmp;
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