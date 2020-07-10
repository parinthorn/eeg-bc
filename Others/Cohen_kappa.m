function k = Cohen_kappa(A1,A2,isNoDiag)
ndiag = 0;
if nargin==3
    if isNoDiag
       ndiag = floor(sqrt(length(A1)));
       if (ndiag^2~=length(A1))
            error('sqrt error')
       end
    else
        error('unknown')
    end
else

end
n11=length(find((A1~=0)&(A2~=0)));
n12=length(find((A1~=0)&(A2==0)));
n21=length(find((A1==0)&(A2~=0)));
n22=length(find((A1==0)&(A2==0)));
p=length(A1);
if (n11==p-ndiag)&&(n22==ndiag)
    k=-1;
    return;
end
if (n11==0)&&(n22==p)
    k=-1;
    return;
end
if (n12==0)&&(n21==0)
    k=1;
    return;
end
PrA = (n11+n22)/p;
PrE = ((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22))/p^2;
k = (PrA-PrE)/(1-PrE);
end