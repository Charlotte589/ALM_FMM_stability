function [T,U,V,W] = multiplication_tensor(M,P,N,R)
    % General MMT that computes the tensor such that for A mxp, B pxn, the
    % product C is equal to:
    % 
    n = min([P,M,N]);
    
    T = zeros(M*P,P*N,M*N);
    
    if exist('R','var')
        U = zeros(M*P,R);
        V = zeros(N*P,R);
        W = zeros(M*N,R);
        r = 1;
        for i=1:M
            for j=1:P
                for k=1:N
                    if r <= R
                        U(i+(j-1)*M,r) = 1;
                        V(j+(k-1)*P,r) = 1;
                        W(k+(i-1)*N,r) = 1;
                        r = r +1;
                    end
                end
            end
        end
    end

    
    for i=1:M
        for j=1:P
            for k=1:N
                %disp([i+(j-1)*M,j+(k-1)*P,k+(i-1)*N])
                T(i+(j-1)*M,j+(k-1)*P,k+(i-1)*N) = 1;
            end
        end
    end
end