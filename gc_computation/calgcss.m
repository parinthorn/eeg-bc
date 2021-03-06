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
%   INPUT
%    	A       =   System's state transition matrix
%       C       =   System's output matrix
%       W       =   Covariance matrix of process noise w(t) 
%     	V       =   Covariance matrix of output noise v(t) 
%      	S       =   Covariance matrix of w(t) and v(t)
%   OUTPUT
%     	FR      =   Granger causality via parameter reduction method
%       Sigma   =   process noise covariance matrix in innovation form
%
%	SPECIAL CASE WHEN V = 0 and Ci = 0
%
%	Originally written by NATTAPORN PLUB-IN
%	Revised by ANAWAT NARTKULPAT, 2020
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
        return
    end
    C(rz,:) = [];
    V(rz,:) = [];
    V(:,rz) = [];
    S(:,rz) = [];
    
    %---------------------------- Full model ------------------------------   
    [P,~,~,~] = idare(A',C',W,V,S); % solve RICCATI for full model
    Sigma = C*P*C' + V;
    diagSigma = diag(Sigma);    % collect Sigma_ii of full model
    
    %---------------------------- Reduced model ---------------------------
    %------------- solve RICCATI for all reduced model --------------------
    [mm,~] = size(C);           % size C after remove zeros
    parfor j=1:mm                  % reduce jth parameter
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
        ind(j) = [];                    
        tmpvar{j} = ind;
        tmpvar2{j} = diag(SigmaR);
    end
    for j=1:mm
        FRr(tmpvar{j},j) = log(tmpvar2{j}./diagSigma(tmpvar{j}));
    end
    FR = zeros(m,m); FR(rnz,rnz) = FRr;
end

