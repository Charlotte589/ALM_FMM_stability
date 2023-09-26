function [f,Fn,hn,gr,Un] = get_params_AL(T,x,beta,y,param,opt)
    [f,Fn,hn,gr,~,Un] = get_params_LM(T,x,beta,y,param,opt);
    
end