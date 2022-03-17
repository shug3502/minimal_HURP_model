th = load('mcmc_output_synthetic_data_v212_lambda_linear_in_mnz.mat');
splitTheta = num2cell(th.theta_store((th.burnin+1):th.niter,:,:), [1 2]); %split A keeping dimension 1 and 2 intact
theta_store_plot = vertcat(splitTheta{:});
theta = median(theta_store_plot,1);  
params.lambda = theta(3); %binding rate
params.D_h = theta(2); %10^(-2); %diffusion const for HURP
params.scale = theta(9); %spatial scale for Ran GTP gradient
params.mu = theta(4); %decay rate
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;
params.L = 3.12; %domain size (um)
params.T = 37; %time duration (s)
params.lambda_gtp=params.lambda/theta(10); %allow preferential binding to gdp tubulin
params.lambda_gdp=params.lambda; % only part to change is the binding in the MNZ/GTP cap region
params.x_0=0; %initial position of chromosomes
params.nx=501; %number of points in spatial discretization
params.gamma1=theta(7);
params.gamma2=theta(8);
params.v_plus=theta(5);
params.v_minus=theta(6);
params.lambda_linear = 1;
font_size = 18;

    x = linspace(0,params.L,params.nx);
    t = linspace(0,params.T,50);
l_h_vec = linspace(0,2.5,51);
hurp_gap=zeros(size(l_h_vec));
for ll = 1:length(l_h_vec)
    ll
    %trailing kinetochore
    params.l_h = l_h_vec(ll);
    params.lambda_gtp=params.lambda/theta(10); % define parameters here
params.lambda_gdp=params.lambda; % only part to change is the binding in the gtp region
params.is_gradient_relative_to_chromosomes=1;
params.gradient_shape = "exponential"; %"flat top", "linear bump", "exponential"
params.v=params.v_plus;
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
params.v = params.v_minus;
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;

fun = @(x,t,u,dudx) pdex1pde(x,t,u,dudx,params); % anonymous function
init_fun_lead = @(x) pdex1ic(x,u(end,:),params);
bc_fun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,params);
%sol = pdepe(m,fun,init_fun_lead,@pdex1bc,x,t);
sol = pdepe(m,fun,init_fun_lead,bc_fun,x,t);
u_lead = sol(:,:,1);
    
%for trailing kinetochore sister
params.lambda_gtp=params.lambda/theta(10); % define parameters here
params.lambda_gdp=params.lambda; % only part to change is the binding in the gtp region
params.v = params.v_plus;
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;


init_fun_trail = @(x) pdex1ic(x,u_lead(end,:),params);
fun = @(x,t,u,dudx) pdex1pde(x,t,u,dudx,params); % anonymous function
bc_fun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,params);
%sol = pdepe(m,fun,init_fun_trail,@pdex1bc,x,t);
sol = pdepe(m,fun,init_fun_trail,bc_fun,x,t);
u_trail = sol(:,:,1);
    [~,I] = max(u_trail(end,:));
    hurp_gap(ll) = x(I);
end

figure;
plot(l_h_vec,hurp_gap,'Linewidth',3,'Color','k');
hold all;
plot(linspace(0,2.5,201),1.56*ones(201,1),'k--','Linewidth',1.5);
xline(quantile(theta_store_plot(:,1),0.025),'-.r','','Linewidth',1.5)
xline(quantile(theta_store_plot(:,1),0.975),'-.r','','Linewidth',1.5)
set(gca,'fontsize',font_size);
xlabel('Length of MNZ and GTP-cap (um)');
ylabel({'Distance between kinetochore';'and HURP peak (um)'});
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 11 11])
print(sprintf('plots/HURP_gap_vs_MNZ_gradient_relative_to_chromosomes_%d_%s',...
    params.is_gradient_relative_to_chromosomes,params.gradient_shape),'-depsc')


%now sensitivity to speed
params.l_h = theta(1);
speed_vec = linspace(0,0.1,51);
hurp_gap=zeros(size(speed_vec));
for vv = 1:length(speed_vec)
    vv
    %trailing kinetochore
    params.v_plus = speed_vec(vv);
    params.v_minus = -speed_vec(vv);
params.lambda_gtp=params.lambda/theta(10); % define parameters here
params.lambda_gdp=params.lambda; % only part to change is the binding in the gtp region
params.is_gradient_relative_to_chromosomes=1;
params.gradient_shape = "exponential"; %"flat top", "linear bump", "exponential"
params.v=params.v_plus;
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
params.v = params.v_minus;
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;

fun = @(x,t,u,dudx) pdex1pde(x,t,u,dudx,params); % anonymous function
init_fun_lead = @(x) pdex1ic(x,u(end,:),params);
bc_fun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,params);
%sol = pdepe(m,fun,init_fun_lead,@pdex1bc,x,t);
sol = pdepe(m,fun,init_fun_lead,bc_fun,x,t);
u_lead = sol(:,:,1);
    
%for trailing kinetochore sister
params.lambda_gtp=params.lambda/theta(10); % define parameters here
params.lambda_gdp=params.lambda; % only part to change is the binding in the gtp region
params.v = params.v_plus;
params.mu_gtp = params.mu;
params.mu_gdp = params.mu;


init_fun_trail = @(x) pdex1ic(x,u_lead(end,:),params);
fun = @(x,t,u,dudx) pdex1pde(x,t,u,dudx,params); % anonymous function
bc_fun = @(xl,ul,xr,ur,t) pdex1bc(xl,ul,xr,ur,t,params);
%sol = pdepe(m,fun,init_fun_trail,@pdex1bc,x,t);
sol = pdepe(m,fun,init_fun_trail,bc_fun,x,t);
u_trail = sol(:,:,1);
    [~,I] = max(u_trail(end,:));
    hurp_gap(vv) = x(I);
end
% 

figure;
plot(speed_vec,hurp_gap,'Linewidth',3,'Color','k');
hold all;
plot(linspace(0,0.1,201),1.56*ones(201,1),'k--','Linewidth',1.5);
xline(quantile(theta_store_plot(:,5),0.025),'-.r','','Linewidth',1.5) %use quantiles from posterior (loaded)
xline(quantile(theta_store_plot(:,5),0.975),'-.r','','Linewidth',1.5)
set(gca,'fontsize',font_size);
xlabel('Speed (um/s)');
%ylabel('Length of HURP gap (um)');
ylabel({'Distance between kinetochore';'and HURP peak (um)'});
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 11 11])
print(sprintf('plots/HURP_gap_vs_speed_%d_%s',...
    params.is_gradient_relative_to_chromosomes,params.gradient_shape),'-depsc')


%%%%%%%%%%%%%%%%%%%%%%%
function rate = lambda(x,l,params)
lambda_gtp = params.lambda_gtp;
lambda_gdp = params.lambda_gdp;
if x<l
    rate=lambda_gtp + params.lambda_linear*(lambda_gdp-lambda_gtp)*x/l;
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
pl = params.gamma1*ul; %0; %ul;
ql = -1; %1; %0
pr = params.gamma2;%ur;
qr = -1; %0;
end
%----------------------------------------------



