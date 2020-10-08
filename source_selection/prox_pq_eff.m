function prox_x = prox_pq_eff(z,v,PARAMETERS)
% This function computes proximal operator 
% of composite norm pq v*(||x||_p,q)^q
% for p=2, and q=1 or q=0.5
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
p = PARAMETERS(1);
q = PARAMETERS(2);
gLen = PARAMETERS(3);
n = length(z);
gNo = n/gLen;
normz_vect = ((sum((reshape(z,[gLen,gNo])).^2,1))).^(0.5)';
if (p==2 && q==0.5)
gTag = (normz_vect)>1.5*v^(2/3);
ngTag = (normz_vect)<=1.5*v^(2/3);
coeff_num = @(zz) (16.*zz.^(1.5).*(cos(pi/3-(1/3)*(acos(v/4*(3./zz).^(3/2))))).^3);
coeff_den = @(zz) (3*sqrt(3)*v+16.*zz.^(1.5).*cos(pi/3-(1/3)*(acos(v/4*(3./zz).^(3/2)))).^3);
tmp_prox(gTag) = (coeff_num(normz_vect(gTag))./coeff_den(normz_vect(gTag)));
tmp_prox(ngTag) = 0;
tmp_prox = kron(tmp_prox',ones(gLen,1));
try prox_x = tmp_prox.*z;
catch
    error('error occurred')
end
elseif (p==2 && q==1)
    gTag = (normz_vect)>v;
    ngTag = (normz_vect)<=v;
    tmp_prox(gTag) = (1-v./(normz_vect(gTag)));
    tmp_prox(ngTag) = 0;
    tmp_prox = kron(tmp_prox',ones(gLen,1));
    prox_x = tmp_prox.*z;
end
end

