function Lhat = augmentL(L,mhat,variance)
[r,m] = size(L);
if mhat-m<0
    warning('The augmented L is shortened')
end
if nargin<3
    if (mhat<m)
        Lhat = L(:,1:mhat);
        return
    end
else
    Lhat = [L variance*randn(r,mhat-m)];
end
end