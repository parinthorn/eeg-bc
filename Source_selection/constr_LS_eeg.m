function [C,flag] = constr_LS_eeg(V,L,W,C_in)
% This program solve constrained least square (flag=0)
%       min_C ||V-LCW||_F
%       subject to C(C_in==0)==0
% If multiple solution exists, minimum norm solution is returned instead (flag=-1)
[m,n] = size(C_in);
C= zeros(m,n);
nz_ind = find(norms(C_in,2,2));
Lnz = L(:,nz_ind);
if length(nz_ind)>= size(L,1)
    C(nz_ind,:) = Lnz'*(Lnz*Lnz'\V/W); %minimum norm solution, when solution is not unique
    flag = -1;
else
    C(nz_ind,:) = (Lnz'*Lnz)\(Lnz'*V*W')/(W*W'); % well-defined least square case, solve normal eq
    flag = 0;
end



end