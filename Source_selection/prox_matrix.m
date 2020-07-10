function prox_C = prox_matrix(C,v,PARAMETERS)
tmp = C';
tmp = tmp(:);
prox_C = prox_pq_eff(tmp,v,PARAMETERS);
prox_C = reshape(prox_C,[PARAMETERS(3),length(tmp)/PARAMETERS(3)])';
end