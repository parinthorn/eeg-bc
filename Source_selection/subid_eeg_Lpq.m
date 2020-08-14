% 
%   General subspace identification 
%   -------------------------------
%   
%   The algorithm 'subid' identifies deterministic, stochastic 
%   as well as combined state space systems from IO data.
%
%           [A,B,C,D,K,R] = subid(y,u,i);
% 
%   Inputs:
%           y: matrix of measured outputs
%           u: matrix of measured inputs 
%              for stochastic systems, u = []
%           i: number of block rows in Hankel matrices
%              (i * #outputs) is the max. order that can be estimated 
%              Typically: i = 2 * (max order)/(#outputs)
%           
%   Outputs:
%           A,B,C,D,K,R: combined state space system
%           
%                  x_{k+1) = A x_k + B u_k + K e_k        
%                    y_k   = C x_k + D u_k + e_k
%                 cov(e_k) = R
%                 
%           For deterministic systems: K = R = []
%           For stochastic systems:    B = D = []
%
%   Optional:
%
%           [A,B,C,D,K,R,AUX,ss] = subid(y,i,n,AUX,W,sil);
%   
%           n:    optional order estimate (default [])
%                 if not given, the user is prompted for the order
%           AUX:  optional auxilary variable to increase speed (default [])
%           W:    optional weighting flag
%                      SV:    Singular values based algorithm 
%                             (default for systems with input u)
%                      CVA:   Canonical variate based algorithm
%                             (default for systems without input u)
%           ss:   column vector with singular values
%           sil:  when equal to 1 no text output is generated
%           
%   Example:
%   
%           [A,B,C,D,K,R,AUX] = subid(y,u,10,2);
%           for k=3:6
%              [A,B,C,D] = subid(y,u,10,k,AUX);
%           end
%           
%   Reference:
%   
%           Subspace Identification for Linear Systems
%           Theory - Implementation - Applications
%           Peter Van Overschee / Bart De Moor
%           Kluwer Academic Publishers, 1996
%           Stochastic algorithm:   Figure 3.13 page 90 (positive)
%           Combined algorithm:     Figure 4.8 page 131 (robust)
%
%   Copyright:
%   
%           Peter Van Overschee, December 1995
%           peter.vanoverschee@esat.kuleuven.ac.be
%
%

function [sys_est,C_out,ss_out,state] = subid_eeg_Lpq(y,L,ind_Ctrue,i,n,sil)
kappa_cond = 0;
IsFine = 1;
warning on
  
if (nargin < 6);sil = 0;end

if ~sil
    disp(' ');
    disp('   Subspace Identification');
    disp('   -----------------------');
end

% Check the arguments
if (nargin < 4);error('subid needs at least four arguments');end
if (nargin < 5);n = [];end

% Turn the data into row vectors and check
[l,ny] = size(y);if (ny < l);y = y';[l,ny] = size(y);end
if (i < 0);error('Number of block rows should be positive');end
if (l < 0);error('Need a non-empty output vector');end
if ((ny-2*i+1) < (2*l*i));error('Not enough data points');end

% Determine the number of columns in the Hankel matrices
j = ny-2*i+1;

% Compute the R factor
Y = blkhank(y/sqrt(j),2*i,j); 	% Output block Hankel
if ~sil
    disp('      Computing ... R factor');
end
R = qr(Y');
R = triu(R)';               % R factor
clear Y
R = R(1:2*i*l,1:2*i*l); 	% Truncate



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              BEGIN ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% **************************************
%               STEP 1 
% **************************************

% Set up some matrices
Rf = R(l*i+1:2*l*i,:); 	% Future outputs
Rp = [R(1:l*i,:)];        % Past outputs

% The oblique projection:
% Computed as 6.1 on page 166
% obl/Ufp = Yf/Ufp * pinv(Wp/Ufp) * (Wp/Ufp)
% The extra projection on Ufp (Uf perpendicular) tends to give better
% numerical conditioning (see algo page 131)
% And it is needed both for CVA as MOESP

% Ob  = (Rf/Rp)*Rp; which is the same as 
Ob = [Rf(:,1:l*i),zeros(l*i,l*i)];

% **************************************
%               STEP 2 
% **************************************

% Compute the SVD
if ~sil
    disp('      Computing ... SVD');
end
% Compute the matrix WOW
WOW = Ob;
W1i = triu(qr(Rf'));
W1i = W1i(1:l*i,1:l*i)';
WOW = W1i\WOW;
[U,S,V] = svd(WOW);
U = W1i*U; 		% CVA
ss_out = diag(S);
clear V S WOW


% **************************************
%               STEP 3 
% **************************************

% Determine the order from the principle angle
if isempty(n)
  ang = real(acos(ss_out))*180/pi;
  figure(gcf);hold off;subplot;
  bar(1:l*i,ang);
  title('Principal Angles');
  ylabel('degrees');
  xlabel('Order');
  
  n = 0;
  while (n < 1) || (n > l*i-1)
    n = input('      System order ? ');
    if isempty(n);n = -1;end
  end
  if ~sil
    disp(' ')
  end
end

% automatically choose order from the last peak ind - 1/10 of the ind
if n == -1
  ang = real(acos(ss_out))*180/pi;
  angdiff = ang;
  angdiff(2:end) = angdiff(2:end) - angdiff(1:end-1);
  
  angdiff(find(angdiff < 1.5)) = 0;
  [~,peak_ind] = findpeaks(angdiff);
  n = peak_ind(end) - floor(peak_ind(end)/10);
end

U1 = U(:,1:n); 				% Determine U1


% **************************************
%               STEP 4 
% **************************************

% Determine gam and gamm
gam  = U1*diag(sqrt(ss_out(1:n)));  % Gamma_i
gamm = gam(1:l*(i-1),:);        % Gamma_i-1 (Page 36)
% The pseudo inverses
gam_inv  = pinv(gam); 			% Pseudo inverse
gamm_inv = pinv(gamm); 			% Pseudo inverse

% Recover the state
% state = gam_inv*R((2*m+l)*i+1:end,1:(2*m+l)*i)*Q(1:(2*m+l)*i,:);
state = [];

% **************************************
%               STEP 5 
% **************************************

% Determine the matrices A and C
% mydisp(sil,['      Computing ... System matrices A,C (Order ',num2str(n),')']); 
% RhsT = [  gam_inv*R((2*m+l)*i+1:2*(m+l)*i,1:(2*m+l)*i),zeros(n,l) ; ...
%     R(m*i+1:2*m*i,1:(2*m+l)*i+l)]; 
% 
% LhsT = [        gamm_inv*R((2*m+l)*i+l+1:2*(m+l)*i,1:(2*m+l)*i+l) ; ...
%     R((2*m+l)*i+1:(2*m+l)*i+l,1:(2*m+l)*i+l)];

%==========================================================================
% Seperately compute A by least square and C by sparse row method

Rhs = [  gam_inv*R(l*i+1:2*l*i,1:l*i),zeros(n,l) ]; 

LhsA = gamm_inv*R(l*i+l+1:2*l*i,1:l*i+l);
A = LhsA/Rhs;

LhsC = R(l*i+1:l*i+l,1:l*i+l);
% C is solved from min ||LhsC-L*C*Rhs||^2 + lambda*sum||Cj|| where H=L*C

% PRECOMPUTATION

if kappa_cond
    % PRECOMPUTATION
    Timepoints = 1000;
    C_out = Solve_Lpq_regression(L,Rhs,LhsC,50,IsFine);
    [C_out.kappa_ncvx,C_out.ind_k_ncvx] = kappa_selection_timeseries(y,0.5,Timepoints,L,C_out.alpha,n,C_out.C_Lpq);
    
    [C_out.kappa_cvx,C_out.ind_k_cvx] = kappa_selection_timeseries(y,1,Timepoints,L,C_out.alpha,n,C_out.C_L21);
else
    [C_out] = Solve_Lpq_regression(L,Rhs,LhsC,50,IsFine);
end

[~,~,nc] = size(C_out.C_Lpq);

% Choose C by the first C such that setdiff(ind_C_Lpq,ind_Ctrue) +
% setdiff(ind_C_Lpq,ind_Ctrue) is minimized

C_norm_Lpq = sqrt(sum(C_out.C_Lpq.^2,2));
ind_C_Lpq = cell(nc,1);
sumlength_Lpq = zeros(nc,1);
for kk = 1:nc
    ind_C_Lpq{kk} = find(C_norm_Lpq(:,:,kk));
    sumlength_Lpq(kk) = length(setdiff(ind_C_Lpq{kk},ind_Ctrue)) + length(setdiff(ind_Ctrue,ind_C_Lpq{kk}));
end
[~,C_out.ind_chosen_Lpq.cheated] = min(sumlength_Lpq);

C_norm_L21 = sqrt(sum(C_out.C_L21.^2,2));
ind_C_L21 = cell(nc,1);
sumlength_L21 = zeros(nc,1);
for kk = 1:nc
    ind_C_L21{kk} = find(C_norm_L21(:,:,kk));
    sumlength_L21(kk) = length(setdiff(ind_C_L21{kk},ind_Ctrue)) + length(setdiff(ind_Ctrue,ind_C_L21{kk}));
end
[~,C_out.ind_chosen_L21.cheated] = min(sumlength_L21);


%==========================================================================
% res_Lpq = [LhsA;LhsC] - [A;H_Lpq]*Rhs; 			% Residuals
% res_Lpq_CLS = [LhsA;LhsC] - [A;H_Lpq_CLS]*Rhs;
% 
% res_L21 = [LhsA;LhsC] - [A;H_L21]*Rhs; 			
% res_L21_CLS = [LhsA;LhsC] - [A;H_L21_CLS]*Rhs;

% **************************************
%               STEP 7 
% **************************************

H_Lpq = zeros(l,n,nc);
H_Lpq_CLS = zeros(l,n,nc);
H_L21 = zeros(l,n,nc);
H_L21_CLS = zeros(l,n,nc);

W_Lpq = zeros(n,n,nc);
W_Lpq_CLS = zeros(n,n,nc);
S_Lpq = zeros(n,l,nc);
S_Lpq_CLS = zeros(n,l,nc);
E_Lpq = zeros(l,l,nc);
E_Lpq_CLS = zeros(l,l,nc);

W_L21 = zeros(n,n,nc);
W_L21_CLS = zeros(n,n,nc);
S_L21 = zeros(n,l,nc);
S_L21_CLS = zeros(n,l,nc);
E_L21 = zeros(l,l,nc);
E_L21_CLS = zeros(l,l,nc);
for i = 1:nc
    
    H_Lpq(:,:,i) = L*C_out.C_Lpq(:,:,i);
    H_Lpq_CLS(:,:,i) = L*C_out.C_Lpq_CLS(:,:,i);
    H_L21(:,:,i) = L*C_out.C_L21(:,:,i);
    H_L21_CLS(:,:,i) = L*C_out.C_L21_CLS(:,:,i);
    
    res_Lpq = [LhsA;LhsC] - [A;H_Lpq(:,:,i)]*Rhs;          
    res_Lpq_CLS = [LhsA;LhsC] - [A;H_Lpq_CLS(:,:,i)]*Rhs;
    res_L21 = [LhsA;LhsC] - [A;H_L21(:,:,i)]*Rhs;          
    res_L21_CLS = [LhsA;LhsC] - [A;H_L21_CLS(:,:,i)]*Rhs;
    

        % Determine QSR from the residuals
%         if ~sil
%             disp(['      Computing ... System matrices G,L0 (Order ',num2str(n),')']);
%         end
        % Determine the residuals
	cov_res = res_Lpq*res_Lpq'; 		% Covariance	
    cov_res_CLS = res_Lpq_CLS*res_Lpq_CLS';
    
    W_Lpq(:,:,i) = cov_res(1:n,1:n);
    S_Lpq(:,:,i) = cov_res(1:n,n+1:n+l);
    E_Lpq(:,:,i) = cov_res(n+1:n+l,n+1:n+l);

    W_Lpq_CLS(:,:,i) = cov_res_CLS(1:n,1:n);
    S_Lpq_CLS(:,:,i) = cov_res_CLS(1:n,n+1:n+l);
    E_Lpq_CLS(:,:,i) = cov_res_CLS(n+1:n+l,n+1:n+l);
    
    cov_res = res_L21*res_L21'; 			
    cov_res_CLS = res_L21_CLS*res_L21_CLS';
    
    W_L21(:,:,i) = cov_res(1:n,1:n);
    S_L21(:,:,i) = cov_res(1:n,n+1:n+l);
    E_L21(:,:,i) = cov_res(n+1:n+l,n+1:n+l);

    W_L21_CLS(:,:,i) = cov_res_CLS(1:n,1:n);
    S_L21_CLS(:,:,i) = cov_res_CLS(1:n,n+1:n+l);
    E_L21_CLS(:,:,i) = cov_res_CLS(n+1:n+l,n+1:n+l);
    
%         sig = dlyap(A,Q);
%         G = A*sig*H_Lpq' + S;
%         L0 = H_Lpq*sig*H_Lpq' + R_cov;
% 
          % Determine K and Ro
%         if ~sil
%             disp('      Computing ... Riccati solution')
%         end
% 
%         [P,~,~] = idare(A',H',[],-L0,-G);
%         Ro = L0 - H*P*H';
%         K = (G - A*P*H')*(Ro\eye(size(Ro)));
%         [K,Ro] = gl2kr(A,G,H_Lpq,L0);

end

sys_est.A = A;
sys_est.H_Lpq = H_Lpq;
sys_est.H_Lpq_CLS = H_Lpq_CLS;
sys_est.H_L21 = H_L21;
sys_est.H_L21_CLS = H_L21_CLS;
  
%   K,K_CLS,C_out,Ro,Ro_CLS,Q,Q_CLS,R,R_CLS,S,S_CLS,ss_out
% C_out.K = K;
% C_out.K_CLS = K_CLS;
% C_out.Ro = Ro;
% C_out.Ro_CLS = Ro_CLS;

C_out.W_Lpq = W_Lpq;
C_out.W_Lpq_CLS = W_Lpq_CLS;
C_out.W_L21 = W_L21;
C_out.W_L21_CLS = W_L21_CLS;

C_out.S_Lpq = S_Lpq;
C_out.S_Lpq_CLS = S_Lpq_CLS;
C_out.S_L21 = S_L21;
C_out.S_L21_CLS = S_L21_CLS;

C_out.E_Lpq = E_Lpq;
C_out.E_Lpq_CLS = E_Lpq_CLS;
C_out.E_L21 = E_L21;
C_out.E_L21_CLS = E_L21_CLS;

end
