% test_stab_DeepMind_345_R47
folder = cd;
load(strcat(folder,'/sols/345/DeepMind_R47.mat'))
 

 %% (5,4,3)

param.M = 5;
param.P = 4;
param.N = 3;
param.R = 47;

TM = multiplication_tensor(param.M,param.P,param.N);
assert(norm(error_CPD(TM,{U,V,W}))==0);

param
[Q1,E1] = pref_stab_asym({U,V,W},param)

%% (4,5,3)

Ut = transpose_factorm(U,param.M,param.P); 
Vt = transpose_factorm(V,param.P,param.N); 
Wt = transpose_factorm(W,param.N,param.M);

param.M = 4;
param.P = 5;
param.N = 3;

TM = multiplication_tensor(param.M,param.P,param.N);
assert(norm(error_CPD(TM,{Ut,Wt,Vt}))==0);

param
[Q2,E2] = pref_stab_asym({Ut,Wt,Vt},param)

%% (5,3,4)

param.M = 5;
param.P = 3;
param.N = 4;

TM = multiplication_tensor(param.M,param.P,param.N);
assert(norm(error_CPD(TM,{Wt,Vt,Ut}))==0);

param
[Q3,E3] = pref_stab_asym({Wt,Vt,Ut},param)

%% (3,4,5)

param.M = 3;
param.P = 4;
param.N = 5;

param
TM = multiplication_tensor(param.M,param.P,param.N);
assert(norm(error_CPD(TM,{Vt,Ut,Wt}))==0);

[Q4,E4] = pref_stab_asym({Vt,Ut,Wt},param)
Q4
%% (4,3,5)

param.M = 4;
param.P = 3;
param.N = 5;

param
TM = multiplication_tensor(param.M,param.P,param.N);
assert(norm(error_CPD(TM,{V,W,U}))==0);

[Q5,E5] = pref_stab_asym({V,W,U},param)

%% (3,5,4)

param.M = 3;
param.P = 5;
param.N = 4;

param
TM = multiplication_tensor(param.M,param.P,param.N);
assert(norm(error_CPD(TM,{W,U,V}))==0);

[Q6,E6] = pref_stab_asym({W,U,V},param)