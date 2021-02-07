function z = normpq_vec(x,p,q,gLen)
n = length(x);
gNo = n/gLen;
tmp = reshape(x,gLen,gNo);
tmp = (sum(abs(tmp).^(p),1)).^(1/p);
z = sum(tmp.^q);
% z = norm(,q)^q;

end
