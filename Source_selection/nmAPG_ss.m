function [C,history] = nmAPG_ss(L,W,V,alpha,qq,STEP_SIZE,IS_LINESEARCH,varargin)
% This program solve problem in the form of 
%           min_C (1/2)||V-LCW ||_F^2 + alpha * \sum \lambda ||C_i^T||_2^q
% INPUT : C0 (optional)
%% Proximal
pp = 2;
% qq = 0.5;
gLen = size(W,1);
%% Algorithm parameters
d = 0.7;
eta = 0.9;
rho = 0.5;
IT_MAX = 200000;
TOL = 1e-6;
ALLPRINT = 0;FREQ = 1000;
% IS_LINESEARCH = 0;
ac = STEP_SIZE;
as = ac;
%%
if isempty(varargin)
    C0 = zeros(size(L,2), gLen); % REMOVE
else
    C0 = varargin{1};
end
obj = @(X) (1/2)*(norm(V-L*X*W,'fro'))^2+alpha*normpq(X,pp,qq);
k=1;
Ckm1 = C0;
Sk = Ckm1;
Ck = Ckm1;
tk = 1;
tkm1 = 0;
uk = obj(Ck);
qk = 1;
if ALLPRINT
    fprintf('%3s\t%10s\t%10s\n','iter', 'objective','step size');
end

% PRECOMPUTATION
LtL = L'*L;
LtVWt = L'*V*W';
WWt = W*W';

while k<IT_MAX
    ta = tic;
    Yk = Ck+(tkm1/tk)*(Sk-Ck)+((tkm1-1)/tk)*(Ck-Ckm1);
%     GradYk = loss_derivative(L,Yk,Z,W);
    tmp1 = LtL*Yk*WWt;
    %Yk-ay*GradYk;
    if (IS_LINESEARCH) && (k>2)
        dY = Yk-Ykm1;dr = LtL*dY*WWt;
        as = norm(dY,'fro')^2/trace(dY'*dr);
        BACKTRACK = 1;
        while BACKTRACK == 1
            Forward_step1 = Yk+as*(LtVWt-tmp1);
            Skp1 = prox_matrix(Forward_step1,as*alpha,[pp,qq,gLen]);
            objS = obj(Skp1);
            if (objS<=obj(Yk)-d*norm(Skp1-Yk,'fro')^2) || (objS<=uk-d*norm(Skp1-Yk,'fro')^2)
                BACKTRACK =0;
            end
            as =as*rho;
            if as<STEP_SIZE
                as= STEP_SIZE;
                BACKTRACK = 0;
            end
        end
        
    else
        Forward_step1 = Yk+as*(LtVWt-tmp1);
        Skp1 = prox_matrix(Forward_step1,as*alpha,[pp,qq,gLen]);
        objS = obj(Skp1);
    end
    
    if objS <= uk-d*(norm(Skp1-Yk,'fro')^2)
        Xkp1=Skp1;
        objval = objS;
    else
        tmp2 = LtL*Ck*WWt;
        if (IS_LINESEARCH) && (k>2)
            dZ = Ck-Ykm1; dr = LtL*dZ*WWt;
            ac = norm(dZ,'fro')^2/trace(dZ'*dr);
            BACKTRACK = 1;
            while BACKTRACK ==1
                Forward_step2 = Ck+ac*(LtVWt-tmp2);
                PGS = prox_matrix(Forward_step2,ac*alpha,[pp,qq,gLen]);
                objPGS = obj(PGS);
                if objPGS<=uk-d*norm(PGS-Sk,'fro')^2
                    BACKTRACK = 0;
                end
                ac = rho*ac;
                if ac<STEP_SIZE
                    ac = STEP_SIZE;
                    BACKTRACK = 0;
                end
            end
        else
            Forward_step2 = Ck+ac*(LtVWt-tmp2);
            PGS = prox_matrix(Forward_step2,ac*alpha,[pp,qq,gLen]);
            objPGS = obj(PGS);
        end
        
%         objPGS = obj(PGS);
        if objS<=objPGS
            Xkp1 = Skp1;
            objval = objS;
        else
            Xkp1 = PGS;
            objval = objPGS;
        end
    end
    
    tkp1 = 0.5*(sqrt(4*tk^2+1)+1);
    qkp1 = eta*qk+1;
    ukp1 = (eta*qk*uk+objval)/qkp1;
    history.tpi(k,1) = toc(ta);
    history.obj(k,1) = objval;
    ERR = norm(Xkp1-Ck,'fro');
    ERRprev = norm(Ck,'fro');
    history.relstep(k,1) = ERR/ERRprev;
    if ((ERR<(TOL*ERRprev)) || (all(Xkp1==0,'all'))) && k>1 %|| ((k>1) && (abs(objval-history.obj(k-1))<objTOL*objval))
        C = (Xkp1);
        FLAG = 0;
        break
    else
        C = (Xkp1);
        FLAG = -1;
    end
    
    
    k = k+1;
    if ((mod(k,FREQ)==0) && ALLPRINT)
        fprintf('%3d\t%10.4f\t%10.4f\n',k,objval,norm(Xkp1-Ck,'fro')/norm(Ckm1,'fro'))
%         toc(ta);
    end
    % CACHING
    tk = tkp1;
    tkm1 = tk;
    Sk = Skp1;
%   vk = vkp1;
    Ckm1 = Ck;
    Ck = Xkp1;
    qk = qkp1;
    uk = ukp1;
    Ykm1 = Yk;
%   tmp1_old = tmp1;
    
end

if FLAG==-1
    error('Max iteration exceeded')
end
end

function z = normpq(C,p,q)
z = norms(C,p,2);
z = sum(abs(z).^(q),1);
end