function [C,flag] = constr_LS_eeg(V,L,W,nz_ind)
% This program solve constrained least square (flag=0)
%       min_C ||V-LCW||_F
%       subject to C(C_in==0)==0
% If Constrained Least square is well-defined, flag = 0
% If Constrained LS is not well-defined, return regularized constrained LS
% solution flag = -1.
[m,n] = size(C_in);
m = size(L,2);
n = size(W,1);
C= zeros(m,n);
Lnz = L(:,nz_ind);
if length(nz_ind)>= size(L,1)
    %     C(nz_ind,:) = Lnz'*(Lnz*Lnz'\V/W); %minimum norm solution, when solution is not unique
    A = Lnz'*Lnz;%experimental features
    invWWt = (W*W')\eye(n);
    B = 0.01*(invWWt);
    F = Lnz'*V*W'*invWWt;
    C(nz_ind,:) = sylvester(A,B,F);
    %     C(nz_ind,:) = (tmp+0.01*eye(size(tmp,1)))\(Lnz'*V*W')/(W*W'); %experimental features
    flag = -1;
else
    A = Lnz'*Lnz;
    C(nz_ind,:) = (A)\(Lnz'*V*W')/(W*W'); % well-defined least square case, solve normal eq
    flag = 0;
    if rcond(A)<eps % machine precision
        A = Lnz'*Lnz;%experimental features
        invWWt = (W*W')\eye(n);
        B = 0.01*(invWWt);
        F = L'*V*W'*invWWt;
        C(nz_ind,:) = sylvester(A,B,F);
        flag = -2;
    end
end



end
