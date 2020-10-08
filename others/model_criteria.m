function [bic,aicc] = model_criteria(V,L,W,C)
% This function return bic, aicc measure from the constrained Least square
%     min (1/2) ||V-LXW||_F ^2
%     s.t. X(C==0)==0
% Input: V, L, W, C
% Output: bic score, aicc score
% Written by Parinthorn Manomaisaowapak, 2020
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
[ny,N] = size(V);
E = V-L*C*W;
S = E*E'/(N); S=(S+S')/2;
SqrtS = chol(S,'lower');
logdet = sum(2*log(diag(SqrtS)));
np = length(find(C));
bic = N*logdet+N*(ny*log(2*pi)+1)+np*log(N);
aicc =  N*logdet+N*(ny*log(2*pi)+1)+np*2+2*np*(np+1)/(N-np-1);
end
