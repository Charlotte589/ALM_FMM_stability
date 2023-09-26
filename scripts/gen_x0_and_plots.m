
param.M = 3;
param.P = 3;
param.N = 3;

param.R = 23;

TM = multiplication_tensor(param.M,param.P,param.N);

ind = find(TM);
[ind1,ind2,ind3] = ind2sub(size(TM),ind);
ind = [ind1,ind2,ind3];

%%
opt.l = -1;
opt.u = 1;
opt.discr = 0;
opt.I = [];
opt.kmax = 10;
opt.beta = 0;

n = param.M*param.P*param.N;

%% mpn triv 1s, asym

k0 = 1;
kmax = 100;
results = cell(1,kmax-k0+1);
for k=k0:kmax
    ind2 = ind;
    rng(k)
    
    U0 = zeros(size(TM,1),param.R);
    V0 = zeros(size(TM,2),param.R);
    W0 = zeros(size(TM,3),param.R);

    for i=1:n-param.R
       i1 = randperm(length(ind2),2);

       U0(ind2(i1(1),1),i) = 1;
       V0(ind2(i1(1),2),i) = 1;
       W0(ind2(i1(1),3),i) = 1;
       
       U0(ind2(i1(2),1),i) = 1;
       V0(ind2(i1(2),2),i) = 1;
       W0(ind2(i1(2),3),i) = 1;

       ind2(i1,:) = [];
    end

    j=1;
    for i=(n-param.R)+1:param.R

       U0(ind2(j,1),i) = 1;
       V0(ind2(j,2),i) = 1;
       W0(ind2(j,3),i) = 1;

       j = j+1;
    end
    
    U0 = U0+10^(-1)*randn(size(U0));
    V0 = V0+10^(-1)*randn(size(V0));
    W0 = W0+10^(-1)*randn(size(W0));

    out = ALM(TM,{U0,V0,W0},param,opt);
    out.x0 = {U0,V0,W0};
    
    results{k-k0+1} = out;
end

save(sprintf('results_%d%d%d_R%d_ul%d_mpntriv1s_10-2randn_rng%d-%d_new.mat',param.M,param.P,param.N,param.R,opt.u,k0,kmax),'results')

%% Starting points close to a certain PD {U,V,W}, LM, discr
k0 = 1;
kmax = 10;
results = cell(1,kmax-k0+1);

opt.discr = 1;
opt.beta = 10^-1;
folder = cd;
load(strcat(folder,'/sols/333/Rank 23/x_Smirnov_139nzs_rankJ_534.mat'))
U = x{1};
V = x{2};
W = x{3};

for k=k0:kmax
    rng(k)
    disp(k)
    
    U0 = U+10^(-2)*randn(size(TM,1),param.R);
    V0 = V+10^(-2)*randn(size(TM,2),param.R);
    W0 = W+10^(-2)*randn(size(TM,3),param.R);

    out = LM(TM,{U0,V0,W0},param,opt);
    
    results{k-k0+1} = out;
end

save(sprintf('results_%d%d%d_R%d_ul%d_PD+_10-2_rand_rng%d-%d_LM_discr.mat',param.M,param.P,param.N,param.R,opt.u,k0,kmax),'results')

%% Starting points close to zero, CS

opt.discr = 0;
opt.u = 1;
opt.l = -1;
k0 = 1;
opt.kmax = 1000;
kmax = 100;
results = cell(1,kmax-k0+1);
param.S = 2;
param.T = 7;
for k=k0:kmax
    
    U0 = 10^(-2)*randn(size(TM,1),param.R);

    out = ALM(TM,U0,param,opt);
    
    results{k-k0+1} = out;
end

save(sprintf('results_%d%d%d_R%d_10-1randn_rng%d-%d_ALM.mat',param.M,param.P,param.N,param.R,k0,kmax),'results')


%% Plot costs
figure
hold on
for i=1:length(results)
    cost = results{i}.cost;

    semilogy(cost,'b-o','linewidth',1.5,'Markersize',3)
   
    set(gca, 'YScale', 'log')
end
xlabel('Iteration $k$', 'interpreter','latex')
ylabel('$f(x_k)$', 'interpreter','latex')
axis tight
box on

%% Compute rank Jacobian

ranks = zeros(length(results),1);

for i=1:length(results)
    if norm(error_CPD(TM,results{i}.Ubest)) < 10^-10
        ranks(i) = rank(deriv_CPD(results{i}.Ubest),10^-9);
        disp(i)
    end
end

