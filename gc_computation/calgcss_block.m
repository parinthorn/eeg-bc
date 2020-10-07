%========================================================================== 
%
%   This function is for computing Granger causality (GC) for state-space
%   model as
%
%       x(t+1)      =   Ax(t) + w(t)                
%       y(t)        =   Cx(t) + v(t)           
%   
%   The noise covariance matrx
%   [ W     S ]     =   [ E(ww')    E(wv')  ]
%   [ S'    V ]         [ E(vw')    E(vv')  ]
%
%   The causalities are determined between multiple variables in blocks. 
%
%--------------------------------------------------------------------------
%
%   PARAMETERS
%    	A       =   System's state transition matrix
%       C       =   System's output matrix
%       W       =   Covariance matrix of process noise w(t) 
%     	V       =   Covariance matrix of output noise v(t) 
%      	S       =   Covariance matrix of w(t) and v(t)
%       block   =   Variable block indices (1-dimensional cell array)
%                   e.g. {[1 2 3],[4 6],[5 7 8 9],[10]} for 10 dimensions
%     	FR      =   Granger causality via parameter reduction method
%
%	SPECIAL CASE WHEN V = 0 and Ci = 0
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
%%
function [FR,Sigma] = calgcss_block(A,C,W,V,S,block)
    
    [m,n] = size(C);

    if nargin < 6
        block = cell(m,1);
        n_block = m;
        for i = 1:m
            block{i} = i;
        end
    else
        n_block = length(block);
    end
    
    if nargin < 5 || isempty(S), S = zeros(n,m); end;
    assert(size(W,1) == n && size(W,2) == n,'W must be square matrix with size (n,n)');
    assert(size(V,1) == m && size(V,2) == m,'V must be square matrix with size (m,m)');
    assert(size(S,1) == n && size(S,2) == m,'S must be a matrix with size (n,m)');
    
    rz = [];
    for i=1:n_block
        if((norm(C(block{i},:))==0)) 
            rz = [rz i];     % get zero rows of F
        end
    end
    rnz = setdiff(1:n_block,rz);
    if (length(rnz)<=1)
        FR = eye(n_block);
        return
    end
    
    FR = zeros(n_block,n_block); 
    
    for i = rnz
        %---------------------------- Full model ------------------------------   
        [P,~,~,~] = idare(A',C',W,V,S);  % solve RICCATI for full model
        Sigma = C(block{i},:)*P*C(block{i},:)' + V(block{i},block{i});
        SigmaChol = chol(Sigma,'lower');
        logdetSigma = sum(2*log(diag(SigmaChol)));

        %---------------------------- Reduced model ---------------------------
        %------------- solve RICCATI for all reduced model --------------------
        for j = rnz                  % reduce jth block parameter
            if j ~= i
                Creduce = C;
                Vreduce = V;
                Sreduce = S;
                Creduce(block{j},:) = [];      % reduce jth row of C
                Vreduce(block{j},:) = [];      % reduce jth row of V
                Vreduce(:,block{j}) = [];      % reduce jth column of V 
                Sreduce(:,block{j}) = [];      % reduce jth column of S

                block_reduce = zeros(size(block{i}));
                for k = 1:length(block{i})
                    temp = sum(block{i}(k) > block{j});
                    block_reduce(k) = block{i}(k)-temp;
                end

                [Preduce,~,~] = idare(A',Creduce',W,Vreduce,Sreduce); 
                SigmaR = Creduce(block_reduce,:)*Preduce*Creduce(block_reduce,:)'+Vreduce(block_reduce,block_reduce);  
                SigmaRChol = chol(SigmaR,'lower');
                logdetSigmaR = sum(2*log(diag(SigmaRChol)));
                FR(i,j) = logdetSigmaR - logdetSigma;
            end
        end
    end
    FR(1:n_block+1:n_block^2) = 0; % Set diagonal entries to be zeros

end

