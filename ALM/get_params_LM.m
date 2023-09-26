function [f,Fn,hn,gr,H,Un] = get_params_LM(T,x,beta,y,param,opt)

    xvec = cell2vec(x);
    
    % AL objective
    f = objective_LM(x,T,y,beta,opt,param); 
    
    % norm error CPD
    F = error_CPD(T,x,param);
    Fn = norm(F);
    
    % Equality constraint
    h = equality(x,opt);
    hn = norm(h);
    
    H = hessian_LM(x,y,beta,opt,param);
    
    % gradient augmented Lagrangian
    gr = gradient_LM(x,T,y,beta,param,opt);
   
    Un = norm(xvec);
end