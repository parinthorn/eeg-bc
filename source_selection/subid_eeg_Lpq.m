%
%   L_{p,q} regularized stochastic subspace identification
%   -------------------------------------------
%
%   The algorithm identifies  stochastic state space systems with the
%   following form;
%
%               z(t+1) = A*z(t) + w(t)
%                 x(t) = C*z(t) + n(t)
%                 y(t) = L*x(t) + v(t)
%                      = L*C*z(t) + e(t); H = L*C, e(t) = L*n(t) + v(t)
% 
%   The framework follow the original stochastic subspace identification, 
%   ie, obtaining Kalman states is the same as using 'CVA' in
%   original 'subid' by Overchee, until the step to estimate (A,C). 
%   We then modify the least-squares formulation to contain Lpq regularization 
%   to promote C to have sparse rows.
%
%   The function is in the following form;
%           [sys_est,C_out,ss_out] = subid_eeg_Lpq(y,L,i,n,sil)
%
%   Inputs:
%           y:          matrix of measured outputs
%           L:          lead field matrix
%           i:          number of block rows in Hankel matrices
%                       (i * #outputs) is the max. order that can be estimated
%                       Typically: i = 2 * (max order)/(#outputs)
%           n:          if > 0, order of state-space model (number of states)
%                       if = -1, automatically choosing the order
%                       if = [] (default), plot principle angles for
%                                          determining the order by the user
%           sil:        when equal to 1 no text output is generated
%
%   Outputs:
%           sys_est:    contain A and H (=L*C)
%           C_out:      contain C from each choice of regularized parameter
%                       and also contain bic scores and the chosen index
%           ss_out:     singular values
%
%
%   Reference:
%
%           Subspace Identification for Linear Systems
%           Theory - Implementation - Applications
%           Peter Van Overschee / Bart De Moor
%           Kluwer Academic Publishers, 1996
%           Stochastic algorithm:   Figure 3.13 page 90 (positive)
%
%	    Manuscript: 
%		Granger Causality Inference in EEG Source Connectivity Analysis: A State-Space Approach
%		Parinthorn Manomaisaowapak, Anawat Nartkulpat, Jitkomut Songsiri
%
%	    	https://www.biorxiv.org/content/10.1101/2020.10.07.329276v1
%
%   Originally written by Peter Van Overschee
%   Modified by Anawat Nartkulpat, Parinthorn Manomaisaowapak, 2020
%
%   Copyright:
%           Anawat Nartkulpat, Parinthorn Manomaisaowapak
%

function [sys_est,C_out,ss_out] = subid_eeg_Lpq(y,L,i,n,sil)


IsFine = 1;
warning on

if (nargin < 6);sil = 0;end

if ~sil
    disp(' ');
    disp('   Subspace Identification');
    disp('   -----------------------');
end

% Check the arguments
if (nargin < 3);error('subid needs at least four arguments');end
if (nargin < 4);n = [];end

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

% **************************************
%               STEP 5
% **************************************

% Determine the matrices A and C
%==========================================================================
% Seperately compute A by least square and C by sparse row method

Rhs = [  gam_inv*R(l*i+1:2*l*i,1:l*i),zeros(n,l) ];

LhsA = gamm_inv*R(l*i+l+1:2*l*i,1:l*i+l);
A = LhsA/Rhs;

LhsC = R(l*i+1:l*i+l,1:l*i+l);
% C is solved from min ||LhsC-L*C*Rhs||^2 + lambda*sum||Cj|| where H=L*C

% PRECOMPUTATION

[C_out] = Solve_Lpq_regression(L,Rhs,LhsC,50,IsFine);


[~,~,nc] = size(C_out.C_Lpq);

%==========================================================================
% **************************************
%               STEP 7
% **************************************

% Calculate covariance matrices
% XXX_Lpq: denotes the estimate from using Lpq regularization (non-convex)
% XXX_L21: denotes the estimate from using L21 regularization (convex)
% XXX_CLS: denotes the estimate from constrained least-squares (with a
% given sparsity)
% (without 'CLS', it means the estimate is regularized solution)

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

    % Determine the residuals
    res_Lpq = [LhsA;LhsC] - [A;H_Lpq(:,:,i)]*Rhs;
    res_Lpq_CLS = [LhsA;LhsC] - [A;H_Lpq_CLS(:,:,i)]*Rhs;
    res_L21 = [LhsA;LhsC] - [A;H_L21(:,:,i)]*Rhs;
    res_L21_CLS = [LhsA;LhsC] - [A;H_L21_CLS(:,:,i)]*Rhs;

    % Determince the covariance	matrices for Lpq formulation
	cov_res = res_Lpq*res_Lpq';
    cov_res_CLS = res_Lpq_CLS*res_Lpq_CLS';

    W_Lpq(:,:,i) = cov_res(1:n,1:n);
    S_Lpq(:,:,i) = cov_res(1:n,n+1:n+l);
    E_Lpq(:,:,i) = cov_res(n+1:n+l,n+1:n+l);

    W_Lpq_CLS(:,:,i) = cov_res_CLS(1:n,1:n);
    S_Lpq_CLS(:,:,i) = cov_res_CLS(1:n,n+1:n+l);
    E_Lpq_CLS(:,:,i) = cov_res_CLS(n+1:n+l,n+1:n+l);

    % Determince the covariance	matrices for L21 formulation
    cov_res = res_L21*res_L21';
    cov_res_CLS = res_L21_CLS*res_L21_CLS';

    W_L21(:,:,i) = cov_res(1:n,1:n);
    S_L21(:,:,i) = cov_res(1:n,n+1:n+l);
    E_L21(:,:,i) = cov_res(n+1:n+l,n+1:n+l);

    W_L21_CLS(:,:,i) = cov_res_CLS(1:n,1:n);
    S_L21_CLS(:,:,i) = cov_res_CLS(1:n,n+1:n+l);
    E_L21_CLS(:,:,i) = cov_res_CLS(n+1:n+l,n+1:n+l);
end

% Set output variables
sys_est.A = A;
sys_est.H_Lpq = H_Lpq;
sys_est.H_Lpq_CLS = H_Lpq_CLS;
sys_est.H_L21 = H_L21;
sys_est.H_L21_CLS = H_L21_CLS;

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
