function [out] = LM(TM,x,param,opts)
    % Implementation of the Levenberg-Marquardt method based on:
    % METHODS FOR NON-LINEAR LEAST SQUARES PROBLEMS, K. Madsen, H.B. Nielsen, O. Tingleff. 
    % Informatics and Mathematical Modelling, Technical University of Denmark, 2nd Edition (April 2004).
    % TM is the (matrix multiplication) tensor of which we want to find a polyadoc decomposition (PD). 
    % x is the starting point.
    % param is a struct containing the parameters of the problem, such as,
    % m,p,n, and R (and possibly S and T for CS decompositions).
    % opts is a struct with the following fields:
    % opt.I, containing the indices of the elements in the PD that we want
    % to be zero by adding a quadratic penalty (QP) term. If you don't want to add a constraint, choose opt.I = [];
    % opt.discr: if equal to 1, a QP term is added to obtain
    % discrete decompositions.
    % opt.beta: constant that is used as regularization parameter for the
    % QP term(s).
    % opt.kmax: maximal number of iterations.
    % out is also a structure containing, e.g., the evolution of the cost
    % function, gradient, and the (approximate) solution that is
    % obtained.

    gradtol = 10^-15;
    steptol = 10^-20;
    nu = 2;

    J = jacobian(x,opts,param);
    F = error_func(TM,x,opts,param);
    g = real(J'*F);
    H = real(J'*J);

    f = 0.5*norm(F)^2;
    
    n = size(H,1);
    found = (norm(g,Inf) < gradtol);
    k = 1;
    
    if isfield(opts,'mu')
        mu = opts.mu;
    else
        mu = max(diag(H));
    end
    
    out = struct;    
    out.cost = f; 
    out.rho = [];
    out.gradn = norm(g);
    out.mu = mu;
    out.normx = norm(cell2vec(x),'inf');
    
    while (found == 0) && (k < opts.kmax)
        
       step = real((H + mu*speye(n))\(-g));

        if norm(step) < steptol*(norm(step)+steptol)
            found = 1;
        else
            xvec = cell2vec(x);
            xnew = xvec + step;       
            xnew = vec2cell(xnew,param);
            fnew = costf(TM,xnew,opts,param);
            
            rho = f - fnew;
            rho = 2*rho/(step'*(mu*step-real(g)));  

            if rho > 0
                
                x = xnew;
                F = error_func(TM,x,opts,param);
                f = 0.5*norm(F)^2;
                J = jacobian(x,opts,param);
                H = (J'*J);  
                g = (J'*F);

                found = (norm(g,Inf) < gradtol);
                mu = mu*max([1/3,1-(2*rho-1)^3]);
                nu = 2;
                k = k + 1;
                              
                out.cost = [out.cost; f];
                out.gradn = [out.gradn; norm(g)];
                out.rho = [out.rho; rho];
                out.mu = [out.mu; mu];
                out.normx = [out.normx; norm(cell2vec(x),'inf')];

                fprintf(('It: %i, f: %i, gn: %i, Un: %i, mu: %i, stepn: %i \n'),k,f, norm(g),norm(cell2vec(x),'inf'),mu,norm(step))
                
            else
                mu = mu*nu;
                nu = nu*1.5;
            end
        end       
        
    end
    
    out.xbest = x;
    out.fbest = f;
    
end