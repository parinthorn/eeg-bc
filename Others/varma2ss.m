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
%========================================================================== 
function [sys] = varma2ss(A,C,type)
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
    sys = ss(As,Bs,Cs,0,1);

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
    sys = ss(As,Bs,Cs,0,1);
end
end
