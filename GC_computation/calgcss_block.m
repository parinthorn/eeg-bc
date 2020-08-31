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
    
%     remove_block = [];
%     for i = rz
%         remove_block = [remove_block, block{i}];
%     end
%     C(remove_block,:) = [];
%     V(remove_block,:) = [];
%     V(:,remove_block) = [];
%     S(:,remove_block) = [];
    
    FR = zeros(n_block,n_block); 
    
    for i = rnz
        %---------------------------- Full model ------------------------------   
        [P,~,K,~] = idare(A',C',W,V,S);  % solve RICCATI for full model
        K = K';
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
    FR(1:n_block+1:n_block^2) = 0;
%     FR(rnz,rnz) = FRr;
    %%           GC from C(A-KC)K condition
%     B(:,:,1) = eye(n);   
%     for k = 1:(n-1)
%         B(:,:,k+1) = B(:,:,k)*(A-K*C);      % compute (A-KC)^k-1
%         FCBKK(:,:,k+1) = C*B(:,:,k)*K;      % CBK condition       
%     end
%     FCBKr = vecnorm(FCBKK,1,3)/m;           % return norm1 of FCBKrr for every layer k
%     FCBK = zeros(m,m); 
%     FCBK(rnz,rnz) = FCBKr;
%     FCBK(logical(eye(size(FCBK,1)))) = 0;   % remove diagonal or FCBK

end

