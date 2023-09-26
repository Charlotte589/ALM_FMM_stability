function Ut = transpose_factorm(U,n1,n2)

    Ut = zeros(size(U));
    
    for i=1:length(Ut)
        UTi = reshape(U(:,i),n1,n2)';
        Ut(:,i) = UTi(:);
        
    end

end