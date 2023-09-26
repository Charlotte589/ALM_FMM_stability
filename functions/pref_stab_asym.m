function [Q,E] = pref_stab_asym(U,param)
    V = U{2};
    W = U{3};
    U = U{1};
    
    alfa = zeros(1,param.R);
    beta = zeros(1,param.R);
    e = zeros(size(W,1),1);
    a = zeros(1,param.R);
    b = zeros(1,param.R);
    q = zeros(size(W,1),1);
    
    for i=1:param.R
        k = find(abs(U(:,i))>10^(-3));
        alfa(i) = length(k);
        k = find(abs(V(:,i))>10^(-3));
        beta(i) = length(k);
        a(i) = sum(abs(U(:,i)));         
        b(i) = sum(abs(V(:,i)));
    end
    
    for i=1:length(q)
        k = find(abs(W(i,:))>10^(-3));
        q(i) = length(k) + max((alfa(k)+beta(k)));
        e(i) = sum((a.*b).*abs(W(i,:)));
    end
    
    E = max(e);
    Q = max(q);

end