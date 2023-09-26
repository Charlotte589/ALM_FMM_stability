function F = error_CPD(T,x,param)
    
    if iscell(x)
        U = x{1};
        V = x{2};
        W = x{3};
    elseif isstruct(x)
        U = x.U;
        V = x.V;
        W = x.W;
    else
        [U,V,W] = gen_cyclic_factorm(x,param);
    end

    F = 0;
    
    for i=1:size(U,2)     
        F = F + kron(kron(W(:,i),V(:,i)),U(:,i));
    end
    
    F = F - T(:);
       
end