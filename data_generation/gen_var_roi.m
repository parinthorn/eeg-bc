function [M] = gen_var_roi(PARAMETER,opts) %No. Groups, Similarity
% This function generates a sparse vector autoregressive model
% y(t) = A1*y(t-1) + A2*y(t-2) + ... + Ap*y(t-p) + u(t)
%
% 'A' represents AR coefficients A1,A2,...,Ap and is stored as a p-dimensional array
%
% The input arguments are
% CHANGE EVERY THING TO struct PARAMETER(n,m,p,)
% 'PARAMETER.n': dimension of AR coefficient matrices
% 'PARAMETER.m': sub-block size (n%m ==0)
% 'PARAMETER.p': AR model order
% 'PARAMETER.density': the fraction of nonzero entries in AR coefficients
% 'PARAMETER.group_density' : ROI source's density
% 'opts.flag'       : 1  generate VAR that have n/m ROIs
%                   : 0  generate VAR with same structure as opts.S
% 'opts.AUGMENT     : 1 augment |Ap_ij| <THRESH to sign(Ap_ij)*MAX_GAIN
%                   : 0 coefficients are unconstrained
% 'opts.intra_GC'   : 1 GC connections inside ROI.
%                   : 0 No GC connections inside ROI.
% OUTPUT
% size, binaryF, First index of each partition
% The indices of nonzero entries are saved in 'ind_nz'.
% if p = 0, 'y' is simply a random variable. In this case, A is the
% covariance matrix of u with sparse inverse.
%
% Written by Parinthorn Manomaisaowapak, Anawat Nartkulpat, 2020
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
%%
n = PARAMETER.n;
p = PARAMETER.p;

[M,S] =genS(opts,PARAMETER);

if sum(M.block_size)~=n
    error('sum of block size not equal to given dimension')
end
%% Static case
if (p==0)
    S = sparse(2*eye(n)+sign(sprandsym(n,PARAMETER.density)));
    [i,j]=find(S);
    S = S+sparse(ceil(max(0,-min(eig(S))))*eye(n));
    A = S\eye(n); % covariance matrix with sparse inverse
    %             R = chol(phi);
    %             y = R'*randn(n,Num); % y reduces to a random variable with covariance 'phi'
    ind_nz = sub2ind([n n],i,j);
    %             figure;plot_spy(ind_nz,n,'image');
    %             title('correct sparsity');
    return;
end

%% Randomize AR coefficients
MAX_EIG = 1;
diag_ind = find(eye(n));
k = length(diag_ind);
diag_ind3D = kron(n^2*(0:p-1)',ones(k,1))+kron(ones(p,1),diag_ind);
off_ind3D = (1-repmat(eye(n),1,1,p));
A = zeros(n,n,p);
ii = 0;
while MAX_EIG
    ii = ii+1;
    if ii>0
        [~,S] = genS(opts,PARAMETER);
    end
    %                 disp(ii)
    for k=1:p
        A(:,:,k) = 0.1*sprandn(S(:,:,1));
    end
    
    poles = -0.7+2*0.7*rand(n,p); % make the poles inside the unit circle
    characeq = zeros(n,p+1);
    for jj=1:n
        characeq(jj,:) = poly(poles(jj,:)); % each row is [1 -a1 -a2 ... -ap]
    end
    aux = -characeq(:,2:end);
    A(diag_ind3D) = aux(:); % replace the diagonal entries with stable polynomial
    %                 size(A)
    %                 A(find(off_ind3D));
    if opts.AUGMENT
        THRESH = 1;
        MAX_GAIN = 1;
        A(((abs(A)<THRESH)&(A~=0))&(off_ind3D))=MAX_GAIN.*sign(A(((abs(A)<THRESH)&(A~=0))&(off_ind3D))); % need revision
    end
    
    
    AA = zeros(n,n*p);
    for k=1:p
        AA(1:n,k*n-(n-1):k*n) = A(:,:,k);
    end
    AA = sparse([AA ; [eye(n*(p-1)) zeros(n*(p-1),n)]]);
    %                 disp(max(abs(eigs(AA))))
    if (max(abs(eigs(AA))) < 1)
        MAX_EIG = 0;
        disp(['Spectral radius of VAR process ',num2str(max(abs(eigs(AA))))])
    end
end

ind_nz = find(A(:,:,1));
M.A = A;
M.ind_nz = ind_nz;
end

function [M,S] = genS(opts,PARAMETER)
n = PARAMETER.n;
n_roi = PARAMETER.n_roi;
density = PARAMETER.density;
group_density = PARAMETER.group_density;


flag = opts.flag;
intra_GC = opts.intra_GC;
if n_roi ==1
    n_roi = n;
    group_density = density;
end
if mod(n,n_roi)==0
    m = n/(n_roi);
elseif mod(n,n_roi-1)==0
    m = n/(n_roi-1)-1;
else
    m = floor(n/(n_roi-1));
end
tmp_size = m*(n_roi-1);
mtilde = n-m*(n_roi-1);


if (flag)
    if intra_GC
        S = ((sprand(tmp_size/m,tmp_size/m,density)+eye(tmp_size/m))~=0);
        Stilde = eye(mtilde);
    else
        S = full((sprand(tmp_size/m,tmp_size/m,density))~=0);
        Stilde = zeros(mtilde,mtilde);
    end
    GC_ROI = S;
    S = kron(S,ones(m,m));
    S(S~=0) = sprand(length(S(S~=0)),1,group_density);
    A1 = (sprand(n_roi-1,1,density))~=0;
    A2 = (sprand(1,n_roi-1,density))~=0;
    M.GC_ROI = [GC_ROI A1;A2 0];
    M.binary_F =(S~=0);
    M.ind_cluster_source = [1:m:m*(n_roi-1) m*(n_roi-1)+1];
    M.block_size = diff([M.ind_cluster_source n+1]);
    A1 =  kron(A1,ones(m,mtilde));
    A2 =  kron(A2,ones(mtilde,m));
    A1(A1~=0) = sprand(length(A1(A1~=0)),1,group_density);
    A2(A2~=0) = sprand(length(A2(A2~=0)),1,group_density);
    S = full([S A1;A2 Stilde]);
else
    S = full(opts.S);
end
end