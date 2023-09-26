function x = vec2cell(xvec,param)
    
    R = param.R;
    N = param.N;
    
    if isfield(param,'S')
        x = reshape(xvec,N^2,R);
    else   
        M = param.M;
        P = param.P;
        n1 = M*P;
        n2 = P*N;
        n3 = M*N;
        
        x = cell(1,3);
        
        xi = xvec(1:n1*R);

        x{1} = reshape(xi,n1,R);

        xi = xvec(n1*R+1:(n1+n2)*R);

        x{2} = reshape(xi,n2,R);

        xi = xvec((n1+n2)*R+1:end);

        x{3} = reshape(xi,n3,R);
    end
end