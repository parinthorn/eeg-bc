%==========================================================================
%
%	This function converts VARMA parameters in the form of
%
%               y(t)	= A1y(t-1) + A2y(t-2) + ... + Any(t-n) + e(t) +
%                         C1e(t-1) + C2e(t-2) + ... + Cme(t-m)
%
%	where m <= n, to state-space model object in the form of
%
%               x(t+1)  = Ax(t) + Be(t)
%               y(t)    = Cx(t)
%
%	where the output is a MATLAB state-space object.
%
%	INPUTS
%           A       =   A(:,:,i) is AR matrix indexed i
%           C       =   C(:,:,j) is MA matrix indexed j
%           type    =   'Hamilton' (default) or 'Harvey' form
%	OUTPUTS
%           sys     =   MATLAB state-space object
%
%   Written by ANAWAT NARTKULPAT, 2020
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
function [sys] = varma2ss(A,C,type,Ts)
if nargin < 3
    type = 'Hamilton';
end
[iA,~,kA] = size(A);
[~,~,kB] = size(C);

% ---------------------- Hamilton form ------------------------------------
if strcmp(type,'Hamilton')
%	Create top companion dynamic matrix As from AR matrices
    As = zeros(iA*kA,iA*kA);
    for n = 1:kA 
        As(1:iA,iA*(n-1)+1:iA*n) = A(:,:,n);
        if n ~= kA
            As(iA*n+1:iA*(n+1),iA*(n-1)+1:iA*n) = eye(iA);
        end
    end
%	Create output matrix Cs from MA matrices
    Cs = [eye(iA),zeros(iA,iA*kA-iA)];
    for m = 1:kB
        Cs(:,iA*m+1:iA*(m+1)) = C(:,:,m);
    end
%	The input matrix is [I 0]'
    Bs = [eye(iA);zeros(iA*kA-iA,iA)];
    sys = ss(As,Bs,Cs,0,Ts);

% ---------------------- Harvey form --------------------------------------
elseif strcmp(type,'Harvey')
%	Create left companion dynamic matrix As from AR matrices
    As = zeros(iA*kA,iA*kA);
    for n = 1:kA 
        As(iA*(n-1)+1:iA*n,1:iA) = A(:,:,n);
        if n ~= kA
            As(iA*(n-1)+1:iA*n,iA*n+1:iA*(n+1)) = eye(iA);
        end
    end
%	Create input matrix Bs from MA matrices
    Bs = [eye(iA);zeros(iA*kA-iA,iA)];
    for m = 1:kB
        Bs(iA*m+1:iA*(m+1),:) = C(:,:,m);
    end
%	The output matrix is [I 0]
    Cs = [eye(iA),zeros(iA,iA*kA-iA)];
    sys = ss(As,Bs,Cs,0,Ts);
end
end

