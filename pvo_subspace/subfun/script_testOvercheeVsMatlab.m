load iddata3

[A,B,C,D,K,R] = subid(z3.y,z3.u,10);
sys_overchee = ss(A,B,C,D,z3.Ts);
n4opt =  n4sidOptions('N4Weight','CVA');
sys_matlab = n4sid(z3,1:5,n4opt);

compare(z3,sys_overchee,sys_matlab);