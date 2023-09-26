function f = costf(TM,x,opts,param)

    F = error_func(TM,x,opts,param);
    
    f = 0.5*norm(F)^2;
end
