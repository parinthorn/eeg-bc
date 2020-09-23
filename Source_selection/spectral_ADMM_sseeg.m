function [x,Px, history] = spectral_ADMM_sseeg(G, b, P,a1, a2,pp,qq, PARAMETER, rho,epscor,Ts,varargin)
% This program solve the problem
%  \min_{x} (1/2)||Gx-b||_{2}^{2} + a1 * \mynorm{Px}{pp}{qq}{p} + a2 * \mynorm{Px}{pp}{qq}{pK}
% The input parameters are
%           PARAMETER : [n,p,K]
%           rho       : initial ADMM parameter
%           epscor    :
%           Ts        :
%           varargin  : initial solution
% The output are
%           x         : solution of the problem in sparse format
%           Px        : A projected solution by matrix P
%           history   : Optimization history which has the data of
%                  .r_norm     : primal residual norm
%                  .s_norm     : dual residual norm
%                  .eps_pri     : primal residual norm convergence threshold
%                  .eps_dual     : dual residual norm convergence threshold
%                  .Lagrange     : value of augmented Lagrangian
%                  .rho     : ADMM penalty in each iteration
%                  .fit     : the sum square loss of final solution
%                  .tpi     : primal residual norm
% implemented as in http://proceedings.mlr.press/v54/xu17a/xu17a.pdf
% always converge if qq=1, pp=2
% set a2 = 0
rhotilde = rho;
PRINT_RESULT = 1;
FREQ_PRINT = 10;
MAXITERS = 10000;
ABSTOL = 1e-7; %-7
RELTOL = 1e-5; %-5
% epscor = 0.001;
% Ts = 10;

n = PARAMETER(1);
p = PARAMETER(2);
K = PARAMETER(3);

t_start = tic;



% store variables
nn = n; nd = (n^2-n)*p*(K-1); np = n;
Gtb = G'*b;

GtG = sparse(G'*G);
PtP = sparse(P'*P);

A = [P;P]; At = A';
% ADMM solver

if PRINT_RESULT
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n','iter', ...
        'r norm','eps pri','s norm','eps dual', 'objective');
end

% initial start
optargin = size(varargin,2);
if optargin == 0,
    x1 = zeros(nn,1);
else
    x1 = varargin{1};
end

z1 = zeros(np,1); z2 = zeros(np,1);
z1hat = z1;
z2hat = z2;
Px1 = P*x1;
x2 = Px1;
x3 = Px1;
x1old = x1;
x2old = x2;
x3old = x3;
z1old = z1;
z2old = z2;
z1hatold = z1hat;
z2hatold = z2hat;
obj = 0.5*norm(G*x1-b)^2 + a1*normpq(x2,pp,qq,p)+a2*normpq(x3,pp,qq,p*K);
L = chol(GtG+2*rhotilde*PtP,'lower');
L = sparse(L); U = L';
% Lg = sparse(chol(GtG,'lower'));
z1hatk0 = z1hat;
z2hatk0 = z2hat;
z1k0 = z1;
z2k0 = z2;
x2k0 = x2;
x3k0 = x3;
x1k0 = x1;
history.CHECK_DIVERGE = 1;
for k=1:MAXITERS,


    % x1-update
    c = Gtb+ (rhotilde*(x2+x3)-(z1+z2));


    x1 = U \ (L \ c);
    % x2-update

    Px1 = x1;
    x2 = prox_pq_eff(Px1+z1/rhotilde,a1/rhotilde,[pp,qq,p]);

    % x3-update

    x3 = prox_pq_eff(Px1+z2/rhotilde,a2/rhotilde,[pp,qq,p*K]);


    z1 = z1 + rhotilde*(Px1-x2);
    z2 = z2 + rhotilde*(Px1-x3);


    z1hat =  z1old + rhotilde*(Px1-x2old);
    z2hat = z2old + rhotilde*(Px1-x3old);
    if  (k>5) && mod(k,Ts)==0 %&& (history.r_norm(k-1) > history.eps_pri(k-1))
        dH = -[(x1old-x1k0);(x1old-x1k0)];
        zzhat = [z1hatold-z1hatk0;z2hatold-z2hatk0];
        zz = [z1old-z1k0;z2old-z2k0];
        dG = -([x2old-x2k0;x3old-x3k0]);
        tmp1 = (dH'*zzhat);
        akMG = tmp1/(norm(dH,2)^2);
        akSD = norm(zzhat,2)^2/tmp1;
        if 2*akMG>akSD
            ak=akMG;
        else
            ak = akSD-akMG/2;
        end
        bkSD = norm(zz,2)^2/(dG'*zz);
        bkMG = (dG'*zz)/norm(dG,2)^2;
        if 2*bkMG>bkSD
            bk=bkMG;
        else
            bk = bkSD-bkMG/2;
        end
        akcor = dH'*zzhat/(norm(dH,2)*norm(zzhat,2));
        bkcor = dG'*zz/(norm(dG,2)*norm(zz,2));
        if (akcor > epscor) &&  (bkcor>epscor)
            rhotilde = sqrt(ak*bk);
            L = chol(GtG+2*rhotilde*PtP,'lower');
            L = sparse(L); U = L';
        elseif (akcor > epscor) &&  (bkcor<=epscor)
            rhotilde = ak;
            L = chol(GtG+2*rhotilde*PtP,'lower');
            L = sparse(L); U = L';
        elseif (akcor <= epscor) &&  (bkcor>epscor)
            rhotilde = bk;
            L = chol(GtG+2*rhotilde*PtP,'lower');
            L = sparse(L); U = L';
        end


        z1hatk0 = z1hatold;
        z2hatk0 = z2hatold;
        x1k0 = x1old;
        x2k0 = x2old;
        x3k0 = x3old;
        z1k0 = z1old;
        z2k0= z2old;
    end





    % stopping criterion

    obj = 0.5*norm(G*x1-b)^2 + a1*normpq(x2,pp,qq,p)+a2*normpq(x3,pp,qq,p*K);
    history.objval(k) = obj;
    history.r_norm(k) = norm([Px1-x2;Px1-x3]);
    history.s_norm(k)  = norm(rhotilde*([x2 - x2old+x3-x3old]));
    history.eps_pri(k) = sqrt(np+np)*ABSTOL + RELTOL*max(norm(A*x1), norm([x2;x3]));
    history.eps_dual(k)= sqrt(nn)*ABSTOL + RELTOL*norm([z1;z2]);
    history.Lagrange(k) = obj + rhotilde/2*norm([Px1-x2;Px1-x3]+[z1/rhotilde;z2/rhotilde])^2;
    history.rho(k) = rhotilde;
    if (PRINT_RESULT && mod(k,FREQ_PRINT) == 0)
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n',k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), ...
            history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
            history.s_norm(k) < history.eps_dual(k) )
        history.fit = 0.5*norm(G*x1-b)^2;
        history.CHECK_DIVERGE = 0;
        break;
    end

    x1km1 = x1old;
    x2km1 = x2old;
    x3km1 = x3old;
    z1km1 = z1old;
    z2km1 = z2old;
    z1hatkm1 = z1hatold;
    z2hatkm1 = z2hatold;
    x1old = x1;
    x2old = x2;
    x3old = x3;
    z1old = z1;
    z2old = z2;
    z1hatold = z1hat;
    z2hatold = z2hat;
end
if k==MAXITERS
    history.flag = -1;
else
    history.flag = 1;
end

if PRINT_RESULT
    t_end = toc(t_start);
    history.tpi = t_end/k;
end
x2((x3==0)) = 0;
x=x2;
Px = x2;
end
