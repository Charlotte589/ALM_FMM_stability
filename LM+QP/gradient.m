function g = gradient(TM,x,opt)

    g = jacobian(x,opt)'*error_func(TM,x,opt);
    
end