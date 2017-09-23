% Compares two exponential survival functions

addpath('/home/grad/richard/Documents/mallick/prophaz');
% Simulate data
rng(33330424);
n = 100;
x = rand(n,1);
ind = x < .5;
y = zeros(n,1);
y(ind) = gamrnd(1,1/10,sum(ind),1);
y(~ind) = gamrnd(1,1/1,sum(~ind),1);
% Do some censoring
% cens = gamrnd(1,1,n,1);
cens = 1000 * ones(n,1);
ind = y < cens;
ds = double(ind); % 1 if not censored, 0 if censored;
y(~ind) = cens(~ind);

X = table(x);
[~,I] = sort(y);
y = y(I);
ds = ds(I);
X = X(I,:);

Y = [y, ds];

parpool(8);
TreePH_MCMCparalleltemp(Y,X,'nmcmc',1000,'burn',1000,...
    'filepath','../output/sim1/','seed',1,'saveall',1,...
    'sprop_a',10,'sprop_omega',10);

load('../output/sim1/mcmc_id1.mat')
plot(output.llike)
[~,I] = max(output.llike);
thetree = output.Trees{I};
Treeplot(thetree)

plot(output.As)
plot(output.Omegas)
plot(output.As,output.Omegas,'o')    


% Get posterior of the survival function...
% Suppose we know omega, a,theta, and the true partition...
omega = thetree.omega;
a = thetree.a;
thetas = [2.5,.25];
prt = double(X{:,1} > .5) + 1; % The final partition
nsamp = 10; % number of Monte Carlo samples at each tt

A = zeros(size(Y,1),1);
for jj=1:size(Y,1)
    tmp = 0;
    for kk=1:length(thetas)
        tmp = tmp + sum(Y(prt==kk,1) >= Y(jj,1)) * thetas(kk);
    end
    A(jj) = tmp;
end




% Estimate the survival function at the vector tt
tt = linspace(.01,2,25);
GAM = zeros(nsamp,length(tt));
for ii=1:length(tt)
    ii
    thegam = zeros(nsamp,1);
    U = zeros(nsamp,1);
    jj = 1;
    while tt(ii) >= Y(jj,1) && jj < size(Y,1);
%         A = 0;
%         A2 = 0;
%         for kk=1:length(thetas)
%             A = A + sum(Y(prt==kk,1) >= Y(jj,1)) * thetas(kk);
%             A2 = A2 + sum(Y(prt==kk,1) >= Y(jj+1,1)) * thetas(kk);
%         end
        if jj == 1
            thediff = Y(jj,1);
        else
            thediff = Y(jj,1) - Y(jj-1,1);
        end
        thegam = thegam + gamrnd(a*omega*thediff,1/(a + A(jj)),nsamp,1);
        if Y(jj,2) % only do for non-censored observations
            U = U + slicesampler(nsamp,a + A(jj),a + A(jj+1),.5,10);
            %U = U + slicesample(.5,nsamp,'pdf',@(x) fu(x,a + A(jj),a + A(jj+1),0));
        end
        jj = jj + 1;
    end
%     for kk=1:length(thetas)
%         A = A + sum(Y(prt==kk,1) >= Y(jj,1)) * thetas(kk);
%     end
    if jj == 1
        thediff = tt(ii);
    else
        thediff = tt(ii) - Y(jj-1,1);
    end
    delta = gamrnd(a*omega*thediff,1/(a + A(jj)),nsamp,1);
    GAM(:,ii) = thegam + U + delta;
end
    
ii = 2;
plot(tt,mean(exp(-GAM*thetas(ii))))
hold on
plot(tt,min(exp(-GAM*thetas(ii))))
plot(tt,max(exp(-GAM*thetas(ii))))
fplot(@(x) exp(-x*1),[0,2])
hold off

ypart = Y(X{:,1} < .5,:);
fplot(@(x) ptheta_t(x,ypart,thetree.a,thetree.omega,thetree.theta_shape,thetree.theta_rate,0),[0,15])
ypart = Y(X{:,1} >= .5,:);
fplot(@(x) ptheta_t(x,ypart,thetree.a,thetree.omega,thetree.theta_shape,thetree.theta_rate,0),[0,5])

thetapart = @(x) ptheta_t(x,ypart,thetree.a,thetree.a0,thetree.b0,thetree.theta_shape,thetree.theta_rate,0);
thing = slicesample(1,10000,'pdf',thetapart);
histogram(thing)







            