function out = initialize_output_ALM(out,Un,L,Fn,grn,hn,beta,w,nu)

    out.cost = L;
    out.CPD = Fn;
    out.grn = grn;
    out.hn = hn;
    out.beta = beta;
    out.w = w;
    out.nu = nu;
    out.Un = Un;
    out.k_LM = [];

end