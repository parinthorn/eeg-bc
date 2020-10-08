function prox_C = prox_matrix(C,v,PARAMETERS)
% this function evaluate proximal operator in a matrix form with input C
%             prox_C = argmin_X ||X||_{p,q} + 1/(2v)*||X-C||_F^2
% Written by PARINTHORN MANOMAISAOWAPAK, 2020
% Input: C, v
%       PARAMETER(1) : p, only p=2 can be used
%       PARAMETER(2) : q, only q=1, q=0.5 can be used
%       PARAMETER(3) : gLen, length of group, gLen = # column C in our problem.
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
tmp = C';
tmp = tmp(:);
prox_C = prox_pq_eff(tmp,v,PARAMETERS);
prox_C = reshape(prox_C,[PARAMETERS(3),length(tmp)/PARAMETERS(3)])';
end
