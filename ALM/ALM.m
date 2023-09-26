function out = ALM(T,x,param,opt)
    % Implementation of the Augmented Lagrangian (AL) method with inequality
    % constraint on the elements of x based on:
    
    % J. Nocedal, Jorge and S. J Wright, Numerical optimization, Springer, New York, 1999.
    
    % T is the (matrix multiplication) tensor of which we want to find a CPD. 
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
    % opt.kmax: maximal number of outer iterations.
    % out is also a structure containing, e.g., the evolution of the cost
    % function, gradient, and the (approximate) solution that is
    % obtained.

    % Initialization Lagrange Multiplier
    if isfield(x,'ybest') 
        y = x.ybest;
        x = x.Ubest;
        y1 = y{1};
        y2 = y{2};
    else
        xvec = cell2vec(x);

        n = length(xvec);
        y1 = zeros(length(opt.I)+opt.discr*length(xvec),1);

        y2 = 10^(-2)*zeros(n,1);
        y = {y1,y2};
    end
     
    mu = (10^(-1)); % (Conn:0.1) 1/beta, <= 1
    beta = 1/mu;
    params_ALM.tau = 0.25; %(Conn: 0.01) %beta/tau, < 1
    params_ALM.gamma_bar = 0.1; % (Conn: 0.1) < 1 gamma_1
    
    alfa = min(mu,params_ALM.gamma_bar); %0.1 % < 1
    
    % Initialization tolerances
    params_ALM.alfa_w = 3; % (Conn:1)
    params_ALM.beta_w = 3; % (Conn:1) 
    params_ALM.w_bar = 1; % <= 1 w_0 (Conn:1)
    
    w = params_ALM.w_bar*(alfa)^(params_ALM.alfa_w); % 0.1
     
    opt.gradtol = w;
    
    params_ALM.alfa_nu = 0.9; % Conn:0.1 % < min(1,alfa_w) 
    params_ALM.beta_nu = 0.9; % Conn:0.9 % < min(1,beta_w) 
    params_ALM.nu_bar = 1; % Conn:1 % <= 1 
    
    % Tolerance constraint
    nu = params_ALM.nu_bar*(alfa)^(params_ALM.alfa_nu); 
    k = 1;
    
    gradtol_opt = 10^(-14);
    nu_opt = 10^(-13);
    opt.gradtol = w;
    
    [L,Fn,hn,gr,Un] = get_params_AL(T,x,beta,y,param,opt);
    out = struct;       
    out = initialize_output_ALM(out,Un,L,Fn,[],[],beta,w,nu);
    
    if hn <= nu_opt && norm(gr) <= gradtol_opt
        disp('Starting point already satisfies optimality conditions')
        out.Ubest = x;
        out.ybest = y;
        return
    end

    while k < opt.kmax

        [x,out_LM,k_LM] = LM_ALM(T,x,y,beta,param,opt);
        
        xvec = cell2vec(x);
        
        [L,Fn,hn,gr,Un] = get_params_AL(T,x,beta,y,param,opt);

        if (hn) <= nu
            if hn <= nu_opt && norm(gr) <= gradtol_opt
                out.Ubest = x;
                out.Lbest = L;
                out.ybest = y;
                out = update_output_ALM(out,Fn,beta,w,nu,k_LM,out_LM);
                
                out.params_ALM = params_ALM;
                disp('Optimality conditions satisfied')
                return
            end
            
            h = equality(x,opt);
            y1 = y1 + beta*h;
            
            y2 = max(y2-beta*(opt.u-xvec),zeros(size(y2))) - max(beta*(opt.l-xvec)-y2,zeros(size(y2)));
            y_new = {y1,y2};

            y = y_new;
            
            alfa = mu;
            nu = max(nu*mu^params_ALM.beta_nu,nu_opt);
            w = max(w*alfa^params_ALM.beta_w,gradtol_opt);

            opt.gradtol = w;
                
        else
            mu = max(params_ALM.tau*mu,10^(-16));
            beta = 1/mu;
            alfa = mu*params_ALM.gamma_bar;
            nu = max(params_ALM.nu_bar*alfa^params_ALM.beta_nu,nu_opt);
            w = max(gradtol_opt,params_ALM.w_bar*alfa^(params_ALM.beta_w));
            opt.gradtol = w;
        end
               
        fprintf(('It: %i, Fn: %i, hn: %i, grn: %i, Un: %i, beta: %i, w: %i, nu: %i \n'),...
            k,Fn,hn,norm(gr),norm(xvec,'inf'),beta,w,nu)
        
        k = k + 1;
        
        out = update_output_ALM(out,Fn,beta,w,nu,k_LM,out_LM);
    
    end
    
    out.Ubest = x;
    out.Lbest = L;
    out.ybest = y;
    out.params_ALM = params_ALM;
    
    disp('Maximal number of iteration reached')
   
end