%==========================================================================
%
%	This function generates a ground truth state-space object with its
%   GC matrix.
%
%   The state-space model is generated from a VAR model under diagonal 
%   filter. This model is in the following form.
%
%          z(t+1)  =   A*z(t) + w(t)
%          x(t)    =   C*z(t) + n(t) (by our generating method, n(t) = 0)
%
%   The state-space of VAR model is then passed under the generated filter.
%   A new state-space model with parameters A,B,C,D is then obtained from
%   Av,Bv,Cv and Af,Bf,Cf,Df
%
%   This state-space model is then used as a ground truth model for testing
%   GC causality. The GC matrix is then calculated.
%   
%
%	INPUTS
%           PARAMETER   =   struct contained the following fields:
%               m_active    :   number of active sources
%               m_inactive  :   number of inactive sources
%               lag         :   VAR model's order
%               fp          :   order of filter's numerator
%               fq          :   order of filter's denominator
%               density     :   density of non-zero entries in GC matrix
%               sigma_ar    :   covariance of the noise e(t)
%           
%	OUTPUTS
%           SYSTEM      =   struct contained the following fields
%               PARAMETER   :   struct contained the following fields. This
%                               part contain all parameters for ground
%                               truth generation
%                   m           ;   number of all sources (active+inactive)
%                   m_active    ;   number of active sources
%                   m_inactive  ;   number of inactive sources
%                   lag         ;   VAR model's order
%                   fp          ;   order of filter's numerator
%                   fq          ;   order of filter's denominator
%                   density     ;   density of non-zero entries in GC matrix
%                   ind_active  ;   indices of non-zero row in C
%               sys_ss      :   MATLAB state-space object
%               F           :   GC matrix of sys_ss
%               sigma_ar    :   covariance of the noise e(t)
%               sigma_w     :   covariance of the noise w(t) = B*e(t)
%
%========================================================================== 
function [SYSTEM] = gengcss_eeg(PARAMETER,binaryGC)
if nargin<2
    binaryGC=0;
    flag = 'fullGC';
end
m_active = PARAMETER.m_active;
m_inactive = PARAMETER.m_inactive;
m = m_active+m_inactive;
lag = PARAMETER.lag;
% fp = PARAMETER.fp;
% fq = PARAMETER.fq;
filter_zero = PARAMETER.filter_zero;
filter_pole = PARAMETER.filter_pole;
density = PARAMETER.density;
n_source_cluster = PARAMETER.n_source_cluster;
group_density = PARAMETER.group_density;

SYSTEM.PARAMETER.m = m;
SYSTEM.PARAMETER.m_active = m_active;
SYSTEM.PARAMETER.m_inactive = m_inactive;
SYSTEM.PARAMETER.lag = lag;
fp = length(filter_zero);
fq = length(filter_pole);
SYSTEM.PARAMETER.filter_zero=filter_zero;
SYSTEM.PARAMETER.filter_pole=filter_pole;
SYSTEM.PARAMETER.density = density;
SYSTEM.PARAMETER.group_density = group_density;
SYSTEM.PARAMETER.n_source_cluster = n_source_cluster;

%   Generate VAR model parameters

VARPARAMETER.n = m_active;
VARPARAMETER.n_roi = n_source_cluster;
VARPARAMETER.p = lag;
VARPARAMETER.density = density;
VARPARAMETER.group_density = group_density;

opts = struct('flag',1,'AUGMENT',1,'intra_GC',0);
VARstruct = gen_VAR_ROI(VARPARAMETER,opts);
ind_cluster_source = VARstruct.ind_cluster_source;
Avar = VARstruct.A;

%   Convert Avar to state-space
sys_var = varma2ss(Avar,zeros(m_active,m_active));

%   Generate an invertible stable diagonal filter and convert to ss model
filter_tf = gen_diagfilter_eeg(filter_zero,filter_pole,m_active);
sys_filter = ss(filter_tf,'minimal');

%   Obtain a series connected state-space model of sys_var and sys_filter
sys_beforeT = series_ss2ss(sys_var,sys_filter);
T = rand(size(sys_beforeT.A));
sys = ss2ss(sys_beforeT,T);
if binaryGC
    flag = 'binaryGC';
    sigma_w = sys.B*PARAMETER.sigma_ar*sys.B';
    sigma_n = zeros(m_active);
    F = (Avar(:,:,1)~=0);
    indexF = find(F);
    indexGC = VARstruct.ind_nz;
    indexDiag = [1:m_active+1:m_active^2]';
    indexF = setdiff(indexF,indexDiag);
    indexGC = setdiff(indexGC,indexDiag);
else
    flag = 'FullGC';
    [sA,~] = size(sys.A);
    A = sys.A; B = sys.B; C = sys.C; D = sys.D;
    sigma_w = B*PARAMETER.sigma_ar*B';
    sigma_n = zeros(m_active);
    F = calgcss(A,C,sigma_w,sigma_n,zeros(sA,m_active));
    indexF = find(abs(F)>0.01);
    indexGC = VARstruct.ind_nz;
    indexDiag = [1:m_active+1:m_active^2]';
    indexGC = setdiff(indexGC,indexDiag);
end
%   Calculated the GC matrix of sys (before)


if ~isempty(union(setdiff(indexF,indexGC),setdiff(indexGC,indexF)))
    error('Zero pattern of F may not be correct');
else
    F(setdiff((1:m_active^2)',indexF)) = 0;
end

if m_inactive > 0
    m_in_perroi = floor(m_inactive/n_source_cluster); 
    res = m_inactive - m_in_perroi*n_source_cluster;
    
    ind_cluster_source_new = zeros(size(ind_cluster_source));
    ind_cluster_source_new(1) = 1;
    [nC,mC] = size(sys.C); [nD,mD] = size(sys.D);
    
    Ai = sys.A;
    Bi = eye(mC);
    
    Ci = zeros(nC+m_inactive,mC);
    Di = zeros(nC+m_inactive,mC);
    
    for i = 1:n_source_cluster-1
        ind_cluster_source_new(i+1) = ind_cluster_source(i+1) + i*m_in_perroi;
        Ci(ind_cluster_source_new(i):ind_cluster_source_new(i+1)-1,:) = ...
            [sys.C(ind_cluster_source(i):ind_cluster_source(i+1)-1,:);zeros(m_in_perroi,mC)];
    end
    Ci(ind_cluster_source_new(n_source_cluster):m,:) = ...
        [sys.C(ind_cluster_source(n_source_cluster):m_active,:);zeros(m_in_perroi+res,m_active*(lag+fq))];% butter
    
    sys_new = ss(Ai,Bi,Ci,Di,1);
    Fnew = zeros(m,m);
    
    ind_F_ij = 1:m_active;
    ind_F_ij_new = zeros(1,m_active);
    for i = 1:n_source_cluster-1
        ind_F_ij_new(ind_cluster_source(i):ind_cluster_source(i+1)-1) = ...
            ind_F_ij(ind_cluster_source(i):ind_cluster_source(i+1)-1)+(i-1)*m_in_perroi;
    end
    ind_F_ij_new(ind_cluster_source(i+1):m_active) = ...
        ind_F_ij(ind_cluster_source(i+1):m_active)+i*m_in_perroi;
    
    Fnew(ind_F_ij_new,ind_F_ij_new) = F;
    
    active_index = ind_F_ij_new;
else
    active_index = 1:m_active;
    ind_cluster_source_new = 0;
    Fnew = F;
end

SYSTEM.PARAMETER.ind_active = active_index;
SYSTEM.ind_cluster_source = ind_cluster_source_new;
SYSTEM.ind_cluster_source_old = ind_cluster_source;
SYSTEM.sys_ss = sys_new;
SYSTEM.F = Fnew;
SYSTEM.Fold = F;
SYSTEM.sys_ss_old = sys;
SYSTEM.sigma_w = sigma_w;
SYSTEM.indexGC = indexGC;
SYSTEM.flag = flag;




