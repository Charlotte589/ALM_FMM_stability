function [J] = jacobian(x,opts,param)
    
    J = deriv_CPD(x,param);

    Dh = deriv_equality(x,opts);
    
    J = [J;Dh];
    
end