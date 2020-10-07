function [C,flag] = constr_LS_eeg(V,L,W,nz_ind)
% This program solve constrained least square
%       min_C ||V-LCW||_F
%       subject to C(nz_ind,:)==0
% If Constrained Least square is well-defined, flag = 0
% If Constrained Least square is not well-defined, return regularized constrained LS
% solution flag = -1.
% Written by PARINTHORN MANOMAISAOWAPAK, JITKOMUT SONGSIRI, 2020
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

m = size(L,2);
n = size(W,1);
C= zeros(m,n);
Lnz = L(:,nz_ind);
WWt = W*W';
tmp = min(abs(eig(WWt)));
if tmp<1e-10
    error('WWt is not invertible')
end
A = Lnz'*Lnz;

if (length(nz_ind)>= size(L,1)) || (rcond(A)<eps) % Lnz is fat so Lnz'Lnz is not invertible, We solve L2 regularization instead
    invWWt = (W*W')\eye(n);
    B = 0.01*(invWWt); % gamma = 0.01, L2 penalty
    F = Lnz'*V*W'*invWWt;
    C(nz_ind,:) = sylvester(A,B,F);
    flag = -1;
else
    C(nz_ind,:) = (A)\(Lnz'*V*W')/(W*W'); % well-defined least square case, solve normal eq
    flag = 0;
end



end
