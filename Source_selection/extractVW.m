function [LhsC,Rhs] = extractVW(y,L,u,i,n,AUXin,W,sil);
% This function only find V, W matrix from subid_eeg_Lpq.m
% only for testing
%% Written by PARINTHORN MANOMAISAOWAPAK
warning on

if (nargin < 8);sil = 0;end

mydisp(sil,' ');
mydisp(sil,'   Subspace Identification');
mydisp(sil,'   -----------------------');

% Check the arguments
if (nargin < 4);error('subid needs at least four arguments');end
if (nargin < 5);n = [];end
if (nargin < 6);AUXin = [];end

% Check if its deterministic or stochastic ID
if isempty(u);   ds_flag = 2; 		% Stochastic
else;           ds_flag = 1; 		% Deterministic
end

% Give W its default value
if (nargin < 7);W = [];end
if isempty(W)
  if (ds_flag == 1); W = 'SV'; 		% Deterministic: default to SV
  else;            W = 'CVA';end 	% Stochastic: default to CVA
end


% Turn the data into row vectors and check
[l,ny] = size(y);if (ny < l);y = y';[l,ny] = size(y);end
if (i < 0);error('Number of block rows should be positive');end
if (l < 0);error('Need a non-empty output vector');end
if (ds_flag == 1)
  [m,nu] = size(u);if (nu < m);u = u';[m,nu] = size(u);end
  if (m < 0);error('Need a non-empty input vector');end
  if (nu ~= ny);error('Number of data points different in input and output');end
else
  m = 0;
end
if ((ny-2*i+1) < (2*l*i));error('Not enough data points');end

% Check the weight to be used
Wn = 0;
if (length(W) == 2)
  if (all(W == 'SV') | all(W == 'sv') | all(W == 'Sv'));
    Wn = 1;
    if (ds_flag == 1);Waux = 2;else;Waux = 3;end
  end
end
if (length(W) == 3)
  if (prod(W == 'CVA') | prod(W == 'cva') | prod(W == 'Cva'));
    Wn = 2;
    if (ds_flag == 1);Waux = 3;else;Waux = 1;end
  end
end
if (Wn == 0);error('W should be SV or CVA');end
W = Wn;

% Determine the number of columns in the Hankel matrices
j = ny-2*i+1;

% Check compatibility of AUXin
if (ds_flag == 1);Uaux = u(1,1);else;Uaux = [];end
[AUXin,Wflag] = chkaux(AUXin,i,Uaux,y(1,1),ds_flag,Waux,sil);


% Compute the R factor
if isempty(AUXin)
  Y = blkhank(y/sqrt(j),2*i,j); 	% Output block Hankel
  mydisp(sil,'      Computing ... R factor');
  if (ds_flag == 1)
    U = blkhank(u/sqrt(j),2*i,j); 	% Input block Hankel
    R = triu(qr([U;Y]'))'; 		% R factor
    clear U Y
  else
    R = triu(qr(Y'))'; 			% R factor
    clear Y
  end
  R = R(1:2*i*(m+l),1:2*i*(m+l)); 	% Truncate
else
  R = AUXin(2:2*i*(m+l)+1,1:2*(m+l)*i);
  bb = 2*i*(m+l)+1;
end


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

mi2  = 2*m*i;
% Set up some matrices
if isempty(AUXin) | (Wflag == 1)
  Rf = R((2*m+l)*i+1:2*(m+l)*i,:); 	% Future outputs
  Rp = [R(1:m*i,:);R(2*m*i+1:(2*m+l)*i,:)]; % Past (inputs and) outputs
  if (ds_flag == 1)
    Ru  = R(m*i+1:2*m*i,1:mi2); 	% Future inputs
    % Perpendicular Future outputs
    Rfp = [Rf(:,1:mi2) - (Rf(:,1:mi2)/Ru)*Ru,Rf(:,mi2+1:2*(m+l)*i)];
    % Perpendicular Past
    Rpp = [Rp(:,1:mi2) - (Rp(:,1:mi2)/Ru)*Ru,Rp(:,mi2+1:2*(m+l)*i)];
  end
end

% The oblique projection:
% Computed as 6.1 on page 166
% obl/Ufp = Yf/Ufp * pinv(Wp/Ufp) * (Wp/Ufp)
% The extra projection on Ufp (Uf perpendicular) tends to give better
% numerical conditioning (see algo page 131)
% And it is needed both for CVA as MOESP

if isempty(AUXin)
  if (ds_flag == 1)
    % Funny rank check (SVD takes too long)
    % This check is needed to avoid rank deficiency warnings
    if (norm(Rpp(:,(2*m+l)*i-2*l:(2*m+l)*i),'fro')) < 1e-10
      Ob  = (Rfp*pinv(Rpp')')*Rp; 	% Oblique projection
    else
      Ob = (Rfp/Rpp)*Rp;
    end
  else
    % Ob  = (Rf/Rp)*Rp; which is the same as
    Ob = [Rf(:,1:l*i),zeros(l*i,l*i)];
  end
else
  % Determine Ob from AUXin
  Ob = AUXin(bb+1:bb+l*i,1:2*(l+m)*i);
  bb = bb+l*i;
end


% **************************************
%               STEP 2
% **************************************

% Compute the SVD
if isempty(AUXin) | (Wflag == 1)
  mydisp(sil,'      Computing ... SVD');
  % Compute the matrix WOW we want to take an SVD of
  % W = 1 (SV), W = 2 (CVA)
  if (ds_flag == 1)
    % Extra projection of Ob on Uf perpendicular
    WOW = [Ob(:,1:mi2) - (Ob(:,1:mi2)/Ru)*Ru,Ob(:,mi2+1:2*(m+l)*i)];
  else
    WOW = Ob;
  end
  if (W == 2)
    W1i = triu(qr(Rf'));
    W1i = W1i(1:l*i,1:l*i)';
    WOW = W1i\WOW;
  end
  [U,S,V] = svd(WOW);
  if W == 2;U = W1i*U;end 		% CVA
  ss = diag(S);
  clear V S WOW
else
  U = AUXin(bb+1:bb+l*i,1:l*i);
  ss = AUXin(bb+1:bb+l*i,l*i+1);
end


% **************************************
%               STEP 3
% **************************************

% Determine the order from the singular values
if isempty(n)
  figure(gcf);hold off;subplot;
  if (W == 2)
    bar([1:l*i],real(acos(ss))*180/pi);
    title('Principal Angles');
    ylabel('degrees');
  else
    H = bar([1:l*i],ss); xx = get(H,'XData'); yy = get(H,'YData');
    semilogy(xx,yy+10^(floor(log10(min(ss)))));
    axis([0,length(ss)+1,10^(floor(log10(min(ss)))),10^(ceil(log10(max(ss))))]);
    title('Singular Values');
  end
  xlabel('Order');
  n = 0;
  while (n < 1) | (n > l*i-1)
    n = input('      System order ? ');
    if isempty(n);n = -1;end
  end
  mydisp(sil,' ')
end

U1 = U(:,1:n); 				% Determine U1


% **************************************
%               STEP 4
% **************************************

% Determine gam and gamm
gam  = U1*diag(sqrt(ss(1:n)));  % Gamma_i
gamm = gam(1:l*(i-1),:);        % Gamma_i-1 (Page 36)
% The pseudo inverses
gam_inv  = pinv(gam); 			% Pseudo inverse
gamm_inv = pinv(gamm); 			% Pseudo inverse


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

Rhs = [  gam_inv*R((2*m+l)*i+1:2*(m+l)*i,1:(2*m+l)*i),zeros(n,l) ; ...
         R(m*i+1:2*m*i,1:(2*m+l)*i+l)];

LhsA = gamm_inv*R((2*m+l)*i+l+1:2*(m+l)*i,1:(2*m+l)*i+l);
A = LhsA/Rhs;

LhsC = R((2*m+l)*i+1:(2*m+l)*i+l,1:(2*m+l)*i+l);
% C is solved from min ||LhsC-L*C*Rhs||^2 + lambda*sum||Cj|| where H=L*C
end
