function xn = normalize_max(x)
    U = x{1};
    V = x{2};
    W = x{3};
    
    Un = U;
    Vn = V;
    Wn = W;
    
    R = size(U,2);
    
    for i=1:R
        n = max(abs(U(:,i)))*max(abs(V(:,i)))*max(abs(W(:,i)));
        Un(:,i) = (n^(1/3))*U(:,i)/max(abs(U(:,i)));
        Vn(:,i) = (n^(1/3))*V(:,i)/max(abs(V(:,i)));
        Wn(:,i) = (n^(1/3))*W(:,i)/max(abs(W(:,i)));
    end
   
    xn = {Un,Vn,Wn};
end