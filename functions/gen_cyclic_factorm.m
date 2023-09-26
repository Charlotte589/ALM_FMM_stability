function [U,V,W] = gen_cyclic_factorm(U,param)

    [A,B,C,D] = subdiv_U_ABCD(U,param);
    
    V = [A D B C];

    W = [A C D B];


end