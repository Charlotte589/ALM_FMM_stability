function out = initialize_output_LM(out,Un,f,Fn,hn,grn)

    out.cost = f;
    out.CPD = Fn;
    out.hn = hn;
    out.grn = grn;
    out.Un = Un;

end