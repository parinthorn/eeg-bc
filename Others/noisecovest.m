
function[N,V,costobj] = noisecovest(E,L,VOPTION)
% 
% USAGE: 
% [N,V,KLobj] = noisecovest(E,L,'homo');
% [N,V,KLobj] = noisecovest(E,L,'hetero');
% 
% original problem: minimize KL(N,V) = (1/2)Tr(inv(E)*(LNL'+V)) - logdet(LNL'+V) + logdet(E)
%         subject to N >= 0 , V >= 0,  Tr(E) = Tr(LNL') + Tr(V)
% 
% We solve the problem when N = n*I and for V, we split in two cases: 
% case 1: VOPTION = 'homo' --> V = v*I
% case 2: VOPTION = 'hetero' --> V = diag(v)
% 
% For 'homo' option, the problem reduces to a scalar problem in variable
% 'n' and we solve using bisection (see math detail in the notebook)
% 
% For 'hetero' option, the variable reduces to vector v only but the convex problem is still an SDP 
% where solving the dual has many variables (if r is large)
% 
% Note: for 'hetero' option, we still solve the convex problem using CVX
% 
% Written by JITKOMUT SONGSIRI

[r,m] = size(L);
LLt = L*L'; LLt = (LLt+LLt)/2;
invE = E\eye(r);
det(invE)

switch VOPTION
    case 'homo'
        a = trace(LLt)/r; b = trace(E)/r; c = trace(L'*invE*L)/2 - a*trace(invE)/2; 
        zerograd = @(x) c-trace(  (x*(LLt-a*eye(r))+ b*eye(r))\(LLt-a*eye(r)) ) ;
%         obj = @(x) c*x - log_det( x*( LLt-a*eye(r)) + b*eye(r) ) + log_det(E) + (b/2)*trace(invE);
        obj = @(x) c*x - log_det( x*( LLt-a*eye(r)) + b*eye(r) ) + (b/2)*trace(invE);

        grad0 = zerograd(0); 
        % when r > m (L is tall) the function at b/a is going to infty and
        % numerical error occur since log_det returns a finite value, so we
        % perturb the point b/a by epsilon
        
        gradba = zerograd(b/a-(1e-7)); 
        
        if (grad0 > 0)&&( gradba > 0)
            n = 0; v = trace(E)/r;
        elseif ( grad0 < 0)&&( gradba < 0)
            n = trace(E)/trace(LLt);
            v = 0;
        elseif sign( grad0) ~= sign( gradba) 
            lb = 0; ub = b/a; % Start Bisection 
            while abs(ub-lb) > 1e-10
                n = (lb+ub)/2;
                if zerograd(n) < 0 
                    lb = n;
                else
                ub = n;
                end
            end
            v = -a*n+b;
        end
        
        N = n*eye(m); V = v*eye(r);
        n
        costobj = obj(n);
        
    case 'hetero'
        
        
        a = 1/trace(LLt); b = trace(E)/trace(LLt); c = (diag(invE)-trace(L'*invE*L)*a)/2;
        
        cvx_begin 
        variable v(r,1)
%         minimize (c'*v-log_det( -a*sum(v)*LLt+diag(v) + b*LLt) + log_det(E) + (b/2)*trace(L'*invE*L))
        minimize (c'*v-log_det( -a*sum(v)*LLt+diag(v) + b*LLt) + (b/2)*trace(L'*invE*L))
        subject to
            v >= 0;
            -a*sum(v)+b >= 0;
        cvx_end
        n = -a*sum(v)+b; 
        
        N = n*eye(m); V = diag(v);
        costobj = cvx_optval;
end

