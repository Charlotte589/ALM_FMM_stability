function E = func_7param_323_R15(param,i)

    a = param(1);
    b = param(2);
    c = param(3);
    d = param(4);
    e = param(5);
    g = param(6);
    h = param(7);
    
    U = zeros(6,15);
    V = U; W = zeros(9,15);
    
    U( 1 , 2 ) =  1 ;
    U( 1 , 5 ) =  1 ;
    U( 1 , 6 ) =  1 ;
    U( 1 , 7 ) =  1 ;
    U( 1 , 11 ) =  1 ;
    U( 1 , 15 ) =  1 ;
    U( 2 , 2 ) =  h ;
    U( 2 , 8 ) =  1 ;
    U( 2 , 10 ) =  1 ;
    U( 2 , 13 ) =  1 ;
    U( 2 , 14 ) =  1 ;
    U( 3 , 4 ) =  1 ;
    U( 3 , 8 ) =  b ;
    U( 3 , 10 ) =  b ;
    U( 3 , 11 ) =  -g ;
    U( 3 , 14 ) =  b ;
    U( 4 , 1 ) =  1 ;
    U( 4 , 3 ) =  1 ;
    U( 4 , 5 ) =  -a ;
    U( 4 , 11 ) =  -a ;
    U( 4 , 15 ) =  -a ;
    U( 5 , 1 ) =  h ;
    U( 5 , 6 ) =  a*h ;
    U( 5 , 8 ) =  -a ;
    U( 5 , 10 ) =  -a ;
    U( 5 , 12 ) =  1 ;
    U( 5 , 13 ) =  -a ;
    U( 6 , 3 ) =  -g ;
    U( 6 , 8 ) =  -a*b ;
    U( 6 , 9 ) =  1 ;
    U( 6 , 11 ) =  a*g ;
    U( 6 , 15 ) =  a*g ;
    
    V( 1 , 3 ) =  d ;
    V( 1 , 4 ) =  1 ;
    V( 1 , 11 ) =  d/a ;
    V( 1 , 14 ) =  -a*c ;
    V( 1 , 15 ) =  d ;
    V( 2 , 3 ) =  d/a ;
    V( 2 , 8 ) =  -c ;
    V( 2 , 9 ) =  1 ;
    V( 2 , 10 ) =  -c ;
    V( 2 , 14 ) =  -c ;
    V( 3 , 2 ) =  -e ;
    V( 3 , 6 ) =  e*a ;
    V( 3 , 10 ) =  a ;
    V( 3 , 12 ) =  a ;
    V( 3 , 13 ) =  1 ;
    V( 3 , 14 ) =  a ;
    V( 4 , 1 ) =  -e ;
    V( 4 , 8 ) =  1 ;
    V( 4 , 10 ) =  1 ;
    V( 4 , 12 ) =  1 ;
    V( 4 , 14 ) =  1 ;
    V( 5 , 2 ) =  1 ;
    V( 5 , 3 ) =  a ;
    V( 5 , 7 ) =  a ;
    V( 5 , 11 ) =  1 ;
    V( 5 , 15 ) =  a ;
    V( 6 , 1 ) =  1 ;
    V( 6 , 3 ) =  1 ;
    V( 6 , 5 ) =  1 ;
    V( 6 , 6 ) =  1 ;
    V( 6 , 7 ) =  1 ;
    V( 6 , 15 ) =  1 ;
    
    W( 1 , 3 ) =  a/d ;
    W( 1 , 7 ) =  -1/d ;
    W( 1 , 9 ) =  g ;
    W( 1 , 15 ) =  1/d ;
    W( 2 , 1 ) =  -1/e ;
    W( 2 , 5 ) =  -1/(e*a) ;
    W( 2 , 6 ) =  1/(e*a) ;
    W( 2 , 12 ) =  -h ;
    W( 3 , 5 ) =  -1/a ;
    W( 3 , 7 ) =  1/a ;
    W( 4 , 4 ) =  -b ;
    W( 4 , 10 ) =  1/(a*c) ;
    W( 4 , 12 ) =  1/c ;
    W( 4 , 14 ) =  -1/(a*c) ;
    W( 5 , 12 ) =  1 ;
    W( 5 , 13 ) =  1 ;
    W( 6 , 2 ) =  1/h ;
    W( 6 , 6 ) =  1/(a*h) ;
    W( 6 , 7 ) =  -1/(a*h) ;
    W( 6 , 13 ) =  e ;
    W( 7 , 4 ) =  1 ;
    W( 7 , 9 ) =  1 ;
    W( 8 , 8 ) =  -1/(a*b) ;
    W( 8 , 9 ) =  c ;
    W( 8 , 10 ) =  1/(a*b) ;
    W( 8 , 13 ) =  -1/b ;
    W( 9 , 4 ) =  -d/a ;
    W( 9 , 5 ) =  -1/(a*g) ;
    W( 9 , 11 ) =  -1/g ;
    W( 9 , 15 ) =  1/(a*g) ;

    %% (3,2,3)
    param2.R = 15;
    param2.M = 3;
    param2.P = 2;
    param2.N = 3;
    
    TM = multiplication_tensor(param2.M,param2.P,param2.N);
    assert(norm(error_CPD(TM,{U,V,W}))<=10^-14);

    param2;
    [Q1,E1] = pref_stab_asym({U,V,W},param2);

    %% (2,3,3)

    param2.M = 2;
    param2.P = 3;
    param2.N = 3;
    
    TM = multiplication_tensor(param2.M,param2.P,param2.N);
    assert(norm(error_CPD(TM,{V,W,U}))<=10^-14);

    param2;
    [Q2,E2] = pref_stab_asym({V,W,U},param2);
    
    %% (3,3,2)

    param2.M = 3;
    param2.P = 3;
    param2.N = 2;

    TM = multiplication_tensor(param2.M,param2.P,param2.N);
    assert(norm(error_CPD(TM,{W,U,V}))<=10^-14);

    param2;
    [Q3,E3] = pref_stab_asym({W,U,V},param2);

    Es = {E1,E2,E3};
    E = Es{i};