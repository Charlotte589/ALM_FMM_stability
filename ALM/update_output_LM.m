function out = update_output_LM(out,Un,f,Fn,hn,grn)

    out.cost = [out.cost;f];
    out.CPD = [out.CPD;Fn];
    out.hn = [out.hn;hn];
    out.grn = [out.grn;grn];
    out.Un = [out.Un;Un];

end