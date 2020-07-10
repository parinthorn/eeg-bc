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
%--------------------------------------------------------------------------
%
%   PARAMETERS
%    	A       =   System's state transition matrix
%       C       =   System's output matrix
%       W       =   Covariance matrix of process noise w(t) 
%     	V       =   Covariance matrix of output noise v(t) 
%      	K       =   Covariance matrix of w(t) and v(t)
%     	FR      =   Granger causality via parameter reduction method
%     	FCBK    =   Granger causality via C(A-KC)K condition [BaS:15]
%
%	SPECIAL CASE WHEN V = 0 and Ci = 0
%
%========================================================================== 
%%
function [FR,Sigma] = calgcss(A,C,W,V,S)
    [m,n] = size(C);
    if nargin < 5 || isempty(S), S = zeros(n,m); end;
    assert(size(W,1) == n && size(W,2) == n,'W must be square matrix with size (n,n)');
    assert(size(V,1) == m && size(V,2) == m,'V must be square matrix with size (m,m)');
    assert(size(S,1) == n && size(S,2) == m,'S must be a matrix with size (n,m)');
    
    rz = [];
    for i=1:m
        if((norm(C(i,:))==0)) 
            rz = [rz i];     % get zero rows of C 
        end                
    end
    rnz = setdiff(1:m,rz);
    if (length(rnz)<=1)
        FR = eye(m);
        FCBK = eye(m);
        return
    end
    C(rz,:) = [];
    V(rz,:) = [];
    V(:,rz) = [];
    S(:,rz) = [];
    
    %---------------------------- Full model ------------------------------   
    [P,~,K,~] = dare(A',C',W,V,S);  % solve RICCATI for full model
%     [P,~,K,~] = idare(A',C',W,V);
    % CHANGE TO IDARE
    K = K';
    Sigma = C*P*C' + V;
    diagSigma = diag(Sigma);    % collect Sigma_ii of full model
    
    %---------------------------- Reduced model ---------------------------
    %------------- solve RICCATI for all reduced model --------------------
    [mm,~] = size(C);           % size C after remove zeros
    for j=1:mm                  % reduce jth parameter
        ind = 1:mm;
        Creduce = C;
        Vreduce = V;
        Sreduce = S;
        Creduce(j,:) = [];      % reduce jth row of C
        Vreduce(j,:) = [];      % reduce jth row of V
        Vreduce(:,j) = [];      % reduce jth column of V 
        Sreduce(:,j) = [];      % reduce jth column of S
        [Preduce,~,~] = dare(A',Creduce',W,Vreduce,Sreduce); 
        SigmaR = Creduce*Preduce*Creduce'+Vreduce;
        diagSigmaR = diag(SigmaR);      % collect Sigma_ii of reduced model
        ind(j) = [];                    
        FRr(ind,j) = log(diagSigmaR./diagSigma(ind));
    end
    FR = zeros(m,m); FR(rnz,rnz) = FRr;

    %%           GC from C(A-KC)K condition
    B(:,:,1) = eye(n);   
    for k = 1:(n-1)
        B(:,:,k+1) = B(:,:,k)*(A-K*C);      % compute (A-KC)^k-1
        FCBKK(:,:,k+1) = C*B(:,:,k)*K;      % CBK condition       
    end
    FCBKr = vecnorm(FCBKK,1,3)/m;           % return norm1 of FCBKrr for every layer k
    FCBK = zeros(m,m); 
    FCBK(rnz,rnz) = FCBKr;
    FCBK(logical(eye(size(FCBK,1)))) = 0;   % remove diagonal or FCBK

end

