%==========================================================================
%
%   This function generate a stable, invertible diagonal filer with
%   its ii element equal to
%                   Gii(z)  = pii(z)/qii(z)
%
%                              (z-z_1)(z-z_2)(...)(z-z_2p)
%                           =  ---------------------------
%                              (z-p_1)(z-p_2)(...)(z-p_nq)
%
%   where the order of the polynomial pii(z) and qii(z) are np and nq
%   respectively.
%
%   INPUTS
%           np      =	order of the numerator polynomial
%           nq      =   order of the denominator polynomial
%           N       =   filter matrix size
%           option  =   'same' for using the same poles on all entries
%                       'diff' for random a new set of poles for each entry
%                       no input will uses 'same' as the default
%   OUTPUTS
%           Gz      =   filter's transfer function matrix
%
% Written by Anawat Nartkulpat, 2020
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
%==========================================================================
function Gz = gen_diagfilter_eeg(Z,P,N)
Ts = 1/100;
num = repmat({[0]},N,N); den = repmat({[1]},N,N);
Gz = tf(num,den,Ts);
for i = 1:N
    ZP = zpk(Z,P,1,Ts);
    Gz(i,i) = tf(ZP);
end
end

