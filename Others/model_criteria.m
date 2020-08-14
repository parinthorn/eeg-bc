function [bic,aicc,C,flag] = model_criteria(V,L,ind_nz,W)
% This function solve the problem
%           \argmin_{C} (1/2) ||V-L*C*W||_F^2
%              subject to C(C_in==0)== 0
% If there exists many solution, minimum norm solution is provided.
[C,flag] = constr_LS_eeg(V,L,W,ind_nz);

[ny,N] = size(V);
E = V-L*C*W;
S = E*E'/(N); S=(S+S')/2;
SqrtS = chol(S,'lower');
logdet = sum(2*log(diag(SqrtS)));
np = length(find(C));
bic = N*logdet+N*(ny*log(2*pi)+1)+np*log(N);
aicc =  N*logdet+N*(ny*log(2*pi)+1)+np*2+2*np*(np+1)/(N-np-1);
end
