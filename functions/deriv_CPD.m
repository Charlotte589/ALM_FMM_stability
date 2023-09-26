function J = deriv_CPD(x,param)

    if isstruct(x) || iscell(x)
        if isstruct(x)
            xs = x;
            U = xs.U;
            V = xs.V;
            W = xs.W;  
        else 
            U = x{1};
            V = x{2};
            W = x{3};
        end
        
        n1 = size(U,1);
        n2 = size(V,1);
        n3 = size(W,1);

        R = size(U,2);

        J = cell(1,3*R);

        I1 = eye(n1);
        I2 = eye(n2);
        I3 = eye(n3);

        for i = 1:R
            J{i} = kron(kron(W(:,i),V(:,i)),I1);
            J{R+i} = kron(kron(W(:,i),I2),U(:,i));
            J{2*R+i} = kron(kron(I3,V(:,i)),U(:,i));
        end
    else
        J = cell(1,param.R);
    
        [A,B,C,D] = subdiv_U_ABCD(x,param);

        I = eye(param.N^2);

        for i = 1:param.S
            J{1,i} = kron(kron(I,A(:,i)),A(:,i)) + ...
                kron(kron(A(:,i),I),A(:,i)) + ...
                kron(kron(A(:,i),A(:,i)),I);
        end
        for i = 1:param.T
            J{1,param.S+i} = kron(kron(I,C(:,i)),D(:,i)) + ...
                kron(kron(D(:,i),I),C(:,i)) + ...
                kron(kron(C(:,i),D(:,i)),I);
            J{1,param.S+param.T+i} = kron(kron(I,D(:,i)),B(:,i)) + ...
                kron(kron(B(:,i),I),D(:,i)) + ...
                kron(kron(D(:,i),B(:,i)),I);
            J{1,param.S+2*param.T+i} = kron(kron(I,B(:,i)),C(:,i)) + ...
                kron(kron(C(:,i),I),B(:,i)) + ...
                kron(kron(B(:,i),C(:,i)),I);
        end
    end
    
    J = cell2mat(J);
end