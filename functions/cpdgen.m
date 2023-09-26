function T = cpdgen(x)

    U = x{1};
    V = x{2};
    W = x{3};
    
    T = zeros(size(U,1),size(V,1),size(W,1));
    
   for i=1:size(U,2)
      T = T + reshape(kron(kron(W(:,i),V(:,i)),U(:,i)),size(T)); 
   end

end
