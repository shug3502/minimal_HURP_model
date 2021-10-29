rng(123);
run_mcmc = true;
niter = 10^3;
identifier = "v132_robin";
burnin=niter/2;
nparams=8;

trueparams.lambda = 0.05; %binding rate
trueparams.D_h = 0.005; %10^(-2); %diffusion const for HURP
trueparams.scale = 1.0; %spatial scale for Ran GTP gradient
trueparams.mu = 0.05; %decay rate
trueparams.L = 3; %domain size (um)
trueparams.T = 37; %time duration (s)
trueparams.x_0=0; %initial position of chromosomes
trueparams.nx=51; %number of points in spatial discretization
trueparams.l_h = 0.5;
trueparams.v = 0.03;
trueparams.sigma = 0.01;
trueparams.gamma1 = 0.1;
trueparams.gamma2 = -0.01;

tic; [u_lead_sim,u_trail_sim] = solve_PDE_lead_trail(trueparams); toc;
%add observation noise
u_lead_sim = u_lead_sim + trueparams.sigma*randn(size(u_lead_sim));
u_trail_sim = u_trail_sim + trueparams.sigma*randn(size(u_trail_sim));

color_mat = [0,0,0
    flipud([255,247,243
    253,224,221
    252,197,192
    250,159,181
    247,104,161
    221,52,151
    174,1,126
    122,1,119
    73,0,106]/256)];

font_size = 24;
figure;
x = linspace(0,trueparams.L,trueparams.nx);
subplot(1,2,1);
box on;
hold all;
for j=0:9
    plot(x,u_trail_sim(1+j*5,:),'Color',color_mat(j+1,:),'Linewidth',3);
end
set(gca,'fontsize',font_size);
xlabel('Distance along k-fibre (um)');
ylabel('HURP (a.u.)');
title('Trailing kinetochore (model)');
subplot(1,2,2);
box on;
hold all;
for j=0:9
    plot(x,u_lead_sim(1+j*5,:),'Color',color_mat(j+1,:),'Linewidth',3);
end
set(gca,'fontsize',font_size);
xlabel('Distance along k-fibre (um)');
ylabel('HURP (a.u.)');
title('Leading kinetochore (model)');
%lgd = legend(string(0:9),'Location','west','Orientation','vertical');

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 42 21])
print(sprintf('HURP_PDE_noisy_simulated_data_%s.eps',identifier),'-depsc');

%%%%%%%%%%%%%%%%%%%%%
%Now setup MCMC to infer parameters of interest
%let theta=[l_h,D_h,lambda,mu]; assume we have good estimates of other
%parameters. Fix these at true values.

if run_mcmc
%    theta0 = [0.5,0.001,0.1,0.1,0.05,0.02];
    theta0 = [0.5,0.001,0.1,0.1,0.05,0.1,-0.01,1.0];
    theta_store = NaN(niter,nparams);
    theta = theta0;
    total_acceptances = 0; total_proposals = 0;
    S = diag(0.01*ones(nparams,1));
    loglik = -Inf;
    
    params=trueparams;
    
    while total_acceptances<niter
        %propose new parameters
        theta_star = proposal(theta,S);
        total_proposals = total_proposals + 1;
        
        %solve PDE
        params.l_h = theta_star(1);
        params.D_h = theta_star(2);
        params.lambda = theta_star(3);
        params.mu = theta_star(4);
        params.v = theta_star(5);
        params.gamma1 = theta_star(6);
	params.gamma2 = theta_star(7);
	params.scale = theta_star(8);
%        params.sigma = theta_star(9);
       [u_lead,u_trail] = solve_PDE_lead_trail(params);
        
        %evaluate likelihood
        mask_t = 1+(0:9)*5; %subset of time points to evaluate likelihood at
        mask_x = 1+(0:10)*((trueparams.nx-1)/10); %subset of positions to evaluate likelihood at
        loglik_star = loglikelihood(u_lead(mask_t,mask_x),u_lead_sim(mask_t,mask_x),params) + ...
            loglikelihood(u_trail(mask_t,mask_x),u_trail_sim(mask_t,mask_x),params);
        acceptance_ratio = (loglik_star - loglik) + (prior(theta_star,nparams) - prior(theta,nparams));
        if log(rand(1)) < acceptance_ratio
            %accept
            theta = theta_star;
            loglik=loglik_star;
            total_acceptances = total_acceptances + 1;
            theta_store(total_acceptances,:) = theta;
            if mod(total_acceptances,10)==0
                fprintf(sprintf('iter %d: accept %f - l_h=%f, D_h=%f,lambda=%f, mu=%f, v=%f, gamma1=%f, gamma2=%f, scale=%f, sigma=%f;\n',...
                    total_acceptances,total_acceptances/total_proposals,theta(1),theta(2),theta(3),theta(4),...
		    params.v,params.gamma1,params.gamma2,params.scale,params.sigma));
            end
        end
    end
    fprintf('final acceptance rate: %f \n',total_acceptances/total_proposals);
    fprintf('posterior medians: %f %f %f %f %f %f %f %f\n', median(theta_store((burnin+1):niter,:),1));   
    save(sprintf('mcmc_output_synthetic_data_%s.mat',identifier))
else
    load(sprintf('mcmc_output_synthetic_data_%s.mat',identifier))
end
param_names = {'l_h','D_h','lambda','mu','v','gamma1','gamma2','scale','sigma'};
close all;
figure;
trueparams_vec = [trueparams.l_h,trueparams.D_h,trueparams.lambda,...
    trueparams.mu,trueparams.v,trueparams.gamma1,trueparams.gamma2,trueparams.scale,trueparams.sigma];
for i=1:nparams
    subplot(2,ceil(nparams/2),i);
    histogram(theta_store(:,i),'DisplayStyle','stairs',...
        'Normalization','pdf','EdgeColor','k','LineWidth',2);
    xlabel(param_names{i}); ylabel('Density');
    hold all;
    plot(trueparams_vec(i)*ones(101,1),linspace(0,max(ylim()),101),...
        'r--','linewidth',2);
    set(gca,'fontsize',font_size);
end
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 42 21])
print(sprintf('posterior_histograms_synthetic_data_%s.eps',identifier),'-depsc')
figure;
for i=1:nparams
    subplot(2,ceil(nparams/2),i);
    plot(1:niter,theta_store(:,i),'k','LineWidth',2);
    xlabel('MCMC iteration'); ylabel(param_names{i});
    hold all;
    plot(1:niter,trueparams_vec(i)*ones(niter,1),'r--','linewidth',2);
    set(gca,'fontsize',font_size);
end
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 42 21])
print(sprintf('posterior_traceplot_synthetic_data_%s.eps',identifier),'-depsc')


function theta_star = proposal(theta,S)
%random walk proposal on log space
theta_star = sign(theta).*exp(log(abs(theta)) + mvnrnd(zeros(length(theta),1),S));
end

function p = prior(theta,nparams)
%l_h ~ N(0,1) T[0,];
if theta(1)>=0
    p = log(2*normpdf(theta(1))); %standard normal
else 
    p = -Inf;
end
%D_h ~ N(0,0.1) T[0,];
if theta(2)>=0
    p = p + log(2*normpdf(theta(2),0,0.1));
else 
    p = -Inf;
end
%lambda ~ N(0,1) T[0,];
if theta(3)>=0
    p = p + log(2*normpdf(theta(3),0,1));
else 
    p = -Inf;
end
%mu ~ N(0,1) T[0,];
if theta(4)>=0
    p = p + log(2*normpdf(theta(4),0,1));
else 
    p = -Inf;
end
%v ~ N(0,0.1) T[0,];
if (nparams>=5)
    if (theta(5)>=0)
        p = p + log(2*normpdf(theta(5),0,0.1));
    else 
        p = -Inf;
    end
end
%gamma1 ~ N(0,0.1) T[0,];
if (nparams>=6)
    if (theta(6)>=0)
        p = p + log(2*normpdf(theta(6),0,0.1));
    else 
	p = -Inf;
    end
end
%gamma2 ~ N(0,0.1) T[,0];
if (nparams>=7)
    if (theta(7)<=0)
        p = p + log(2*normpdf(theta(7),0,0.1));
    else 
	p = -Inf;
    end
end
%scale ~ N(0,1) T[0,];
if (nparams>=8)
    if (theta(8)>=0)
        p = p + log(2*normpdf(theta(8),0,1));
    else 
	p = -Inf;
    end
end
%sigma ~ inverse_gamma(a,b)
a = 3; b=0.5;
if (nparams>=9)
    if (theta(9)>=0)
        p = p + log(inversegammapdf(theta(9),a,b));
    else 
        p = -Inf;
    end
end
end
function [ Y ] = inversegampdf( X,A,B )
%inversegampdf Inverse gamma probability density function.
%   Y = inversegampdf(X,A,B) returns the inverse gamma probability density
%   function with shape and scale parameters A and B, respectively, at the
%   values in X. The size of Y is the common size of the input arguments. A
%   scalar input functions is a constant matrix of the same size as the
%   other inputs.
%from https://csdspnest.blogspot.com/2014/03/compute-inverse-gamma-pdf-and-cdf-in.html

Y = B^A/gamma(A)*X.^(-A-1).*exp(-B./X);

end


function [u_lead,u_trail] = solve_PDE_lead_trail(params)
x = linspace(0,params.L,params.nx);
t = linspace(0,params.T,50);
v=params.v; %speed of chromosome movements

%trailing kinetochore
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;
params.lambda_gtp=0; %allow preferential binding to gdp tubulin
params.lambda_gdp=params.lambda; % only part to change is the binding in the MNZ/GTP cap region
params.is_gradient_relative_to_chromosomes=0;
params.gradient_shape = "exponential"; %"flat top", "linear bump", "exponential"
params.v=v; %speed of chromosome movements
m = 0; %symmetry of coordinate system
init_fun = @(x) pdex1ic(x,ones(1,params.nx),params);
fun = @(x,t,u,dudx) pdex1pde(x,t,u,dudx,params); % anonymous function
bc_fun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,params);
%sol = pdepe(m,fun,init_fun,@pdex1bc,x,t);
sol = pdepe(m,fun,init_fun,bc_fun,x,t);
u = sol(:,:,1);

%%%%%%%%%%%%%%%%%%%
%for leading kinetochore sister
params.lambda_gtp=params.lambda;
params.lambda_gdp=params.lambda;
params.v = -v;
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;

fun = @(x,t,u,dudx) pdex1pde(x,t,u,dudx,params);
init_fun_lead = @(x) pdex1ic(x,u(end,:),params);
bc_fun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,params);
%sol = pdepe(m,fun,init_fun_lead,@pdex1bc,x,t);
sol = pdepe(m,fun,init_fun_lead,bc_fun,x,t);
u_lead = sol(:,:,1);

%for trailing kinetochore sister
params.lambda_gtp=0; % define parameters here
params.lambda_gdp=params.lambda; % only part to change is the binding in the gtp region
params.v = v;
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;

init_fun_trail = @(x) pdex1ic(x,u_lead(end,:),params);
fun = @(x,t,u,dudx) pdex1pde(x,t,u,dudx,params);
bc_fun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,params);
%sol = pdepe(m,fun,init_fun_trail,@pdex1bc,x,t);
sol = pdepe(m,fun,init_fun_trail,bc_fun,x,t);
u_trail = sol(:,:,1);
end


%%%%%%%%%%%%%%%%%%%%%%%
function rate = lambda(x,l,params)
lambda_gtp = params.lambda_gtp;
lambda_gdp = params.lambda_gdp;
if x<l
    rate=lambda_gtp;
else
    rate=lambda_gdp;
end
end
function rate = decay_rate(x,l,params)
mu_gtp = params.mu_gtp;
mu_gdp = params.mu_gdp;
if x<l
    rate=mu_gtp;
else
    rate=mu_gdp;
end
end
function f = func(t)
if (length(t)==1) && (t>0)
    f = exp(-1./t);
elseif (length(t)==1) && (t<=0)
    f = 0;
else
    f = (t>0).*exp(-1./t);
end
end
function r = ran_gradient(t,params)
switch params.gradient_shape
    case "flat top"
        %https://math.stackexchange.com/questions/101480/are-there-other-kinds-of-bump-functions-than-e-frac1x2-1
        g = @(t) func(t)./(func(t)+func(1-t));
        h = @(t) g(t-params.scale);
        k = @(t) h(abs(t));
        r = 0.9*(1-k(t))+0.1;
    case "exponential"
        r = exp(-abs(t)/params.scale);
    case "linear bump"
        r = max(2*params.scale-sign(t).*t,0);
    case "kalab fig2c"
        %data extracted from Kalab et al 2006 via the digitize R package
        %linearly interpolate between these points
        dat = [-9.7222222, 0.6904255
            -7.9166667, 0.7063830
            -7.2222222, 0.7063830
            -6.2500000, 0.7117021
            -5.4166667, 0.7202128
            -4.7222222, 0.7308511
            -4.3055556, 0.7329787
            -3.7500000, 0.7351064
            -3.3333333, 0.7446809
            -2.7777778, 0.7585106
            -2.0833333, 0.7723404
            -1.5277778, 0.7882979
            -1.1111111, 0.8074468
            -0.5555556, 0.8255319
            0.0000000, 0.8319149
            0.1388889, 0.8329787
            0.4166667, 0.8308511
            0.8333333, 0.8223404
            0.9722222, 0.8170213
            1.2500000, 0.8095745
            1.5277778, 0.7946809
            2.0833333, 0.7808511
            2.5000000, 0.7648936
            3.0555556, 0.7542553
            3.6111111, 0.7414894
            4.1666667, 0.7340426
            4.8611111, 0.7297872
            5.8333333, 0.7180851
            7.3611111, 0.7063830
            7.6388889, 0.7053191
            8.1944444, 0.7053191
            8.7500000, 0.7000000
            9.7222222, 0.6925532
            10.1388889, 0.6914894];
        r = interp1(dat(:,1),dat(:,2),t);
    case "kalab fig2b"
        dat = [-13.2380952, 225
            -11.5238095, 275
            -9.7142857, 295
            -8.0000000, 325
            -5.9047619, 360
            -4.0952381, 420
            -2.6666667, 460
            -2.0952381, 510
            -1.5238095, 610
            -1.0476190, 720
            -0.5714286, 825
            0.0000000, 870
            0.5714286, 825
            0.7619048, 790
            1.2380952, 685
            1.8095238, 555
            2.7619048, 465
            4.7619048, 415
            7.3333333, 365
            9.6190476, 335
            12.4761905, 275
            15.6190476, 230];
        dat(:,2)=dat(:,2)/max(dat(:,2)); %normalize
        r = interp1(dat(:,1),dat(:,2),t);
end
end
function b = binding_fun(x,u,t,params)
%lambda(x)*g(x) - mu(x)
l_h = params.l_h;
v = params.v*(1-params.is_gradient_relative_to_chromosomes); %to switch between case where gradient moves due to chromosomes moving or not
x_0 = params.x_0;

b = lambda(x,l_h,params)*ran_gradient(x+x_0+v*t,params) - decay_rate(x,l_h,params)*u;
end
function [c,f,s] = pdex1pde(x,t,u,dudx,params) % Equation to solve
D_h = params.D_h; %diffusion const for HURP
v = params.v;
c = 1;
f = D_h*dudx;
s = binding_fun(x,u,t,params) - v*dudx;
end
%----------------------------------------------
function u0 = pdex1ic(xq,u_final,params) % Initial conditions
% if x<0.5
%     u0 = max((0.1*x),0);
% else
%     u0 = max(0.1*(1-x),0);
% end

x = linspace(0,params.L,params.nx);
u0 = interp1(x,u_final,xq);
end
% function u0 = trailic(x,params) % Initial conditions
% hstar = steady_state_of_hurp_model(params);
% u0=hstar(x);
%
% %u0 = exp(-x/0.1);
% end
%----------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t,params) % Boundary conditions
%try robin boundary condition
pl = params.gamma1*ul; %0; %ul;
ql = -1; %1; %0
pr = params.gamma2*ur; %0;%ur;
qr = -1; %1; %0;
end

function L = loglikelihood(data,sim,params)
%likelihood: product over time and space
N = size(data,1)*size(data,2);
L = sum(log(mvnpdf(reshape(data,1,[]),reshape(sim,1,[]),params.sigma*eye(N))),'all');
end
