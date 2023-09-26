function F = error_func(TM,x,opts,param)

    F = error_CPD(TM,x,param);
    
    h = equality(x,opts);
    
    F = [F;h];
    
end