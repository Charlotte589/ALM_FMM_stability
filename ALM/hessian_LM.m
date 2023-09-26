function H = hessian_LM(x,y,beta,opt,param)

    JF = deriv_CPD(x,param);
    xvec = cell2vec(x);
    n = numel(xvec);
    
    Jg = eye(n); % eigenlijk hoort een -1 te staan in de jacobiaan door de bound constraint on the slack variable
    
    [g1,g2] = inequality(x,y{2},beta,opt);
    
    I1 = find(g1);
    I2 = find(g2);
    I = union(I1,I2); % eigenlijk hoort de afgeleide van [l-(g(x)-y_1/beta)]_+ -J_g te zijn (met J_g de jacobiaan van g(x)) 
    % aangezien g(x) in ons geval gewoon gelijk is aan x (U) is de
    % jacobiaan de identiteitsmatrix. De nul rijen worden er met deze find
    % uit gehaald. In de berekening van de hessiaan zouden de -1'en geen
    % verschil maken want door J'*J te nemen worden ze opgehoffen. En in de
    % berekening van de gradient wordt deze matrix niet gebruikt. 
    

    Jh = deriv_equality(x,opt);
    
    J = [JF/sqrt(beta) ; Jg(I,:); Jh];
    
    H = beta*transpose(J)*J;

end