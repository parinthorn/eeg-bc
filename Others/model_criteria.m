function [bic,aicc] = model_criteria(V,L,W,C)
% This function return bic, aicc measure from the constrained Least square
%     min (1/2) ||V-LCW||_F ^2
%     s.t. C(C_Lpq==0)==0
% Input: V, L, W, C
%Output: bic, aicc
% Written by Parinthorn Manomaisaowapak

[ny,N] = size(V);
E = V-L*C*W;
S = E*E'/(N); S=(S+S')/2;
SqrtS = chol(S,'lower');
logdet = sum(2*log(diag(SqrtS)));
np = length(find(C));
bic = N*logdet+N*(ny*log(2*pi)+1)+np*log(N);
aicc =  N*logdet+N*(ny*log(2*pi)+1)+np*2+2*np*(np+1)/(N-np-1);
end
