function xvec = cell2vec(x)

    xvec = [];
    
    if iscell(x)
        for i=1:length(x)
            xi = x{i};
            xvec = [xvec; xi(:)];      
        end
    else
        xvec = x(:);
    end
    
end