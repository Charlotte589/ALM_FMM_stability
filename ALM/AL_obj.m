function L = AL_obj(T,x,y,beta,opt,param)
    % h and y are cells of S + T elements (RxR matrices)
    y1 = y{1};
    y2 = y{2};
    
    f = objective_LM(x,T,y,beta,opt,param);
    
%     L = f - (1/(2*beta))*norm([y1;y2],2)^2;
    L = f;
   
end