%% Solve income fluctuation problem 
% Author: Alessandro Di Nola
% Purpose: solve the income fluctuation problem using golden section search 
% and fast interpolation.
% I wrote my own interpolation routine that is *much faster* than Matlab's
% built-in interp1 or griddedinterpolant.
% See the website quantecon for details about the model and the calibration. 
% -------------------------------------------------------------------------
% Copyright © 2018 by Alessandro Di Nola. All rights 
% reserved.
% -------------------------------------------------------------------------

clear;clc;close all

%% Economic parameters

sigma    = 2;    % CRRA
r        = 0.01; % interest rate (net)
beta     = 0.96; % time discount factor
b        = 0;    % lower bound for asset holdings a
grid_max = 16;   % Upper bound for asset holdings a
R        = 1+r;  % Gross interest rate

% Stochastic process for income shock:
PZ       = [0.60 0.40;  % Markov chain income shock
            0.05 0.95];
z_vals   = [0.5 1.0];   % Income shock values
nz       = length(z_vals);

% Assets Grid
na         = 500;
asset_grid = linspace(-b,grid_max,na)';

% Utility function
if sigma==1
    u = @(c) log(c);
else
    u  = @(c) c^(1-sigma)/(1-sigma);
end

%% Computational parameters

do_plots    = 1;    % flag 0/1
verbose     = 0;    % verbosity level
tiny        = 1e-8; % lower bound on consumption
tol         = 1e-6; % tolerance for VFI
max_iter    = 500;  % maximum num. of iterations 
par.tol_gss = 1e-5; % tolerance for GSS
damp        = 0.0;  % dampening - weight given to previous iterate 
equi_grid   = 1;    % 0=non-equispaced grid, 1=equispaced grid

%% Pack some inputs
par.R          = R;
par.PZ         = PZ;
par.beta       = beta;
par.b          = b;
par.asset_grid = asset_grid;
par.z_vals     = z_vals;
par.na         = na;
par.nz         = nz;
par.tiny       = tiny;
par.u          = u;
par.equi_grid  = equi_grid;

%% Initialization
%Initial conditions for value function and policy function

[v0,c0] = initialize(par);

%% Value function iteration

v_star = compute_fixed_point(v0,tol,max_iter,damp,'bellman',verbose,par);
fprintf(' \n');

%Compute c_pol using final value function
[~,c_pol] = bellman_operator(v_star,par);

%Compute a_pol using c_pol
a_pol = zeros(na,nz);
for iz=1:nz
   a_pol(:,iz) = R*asset_grid+z_vals(iz)-c_pol(:,iz); 
end

%% Some plots

if do_plots==1
% Plot value function
figure(1)
plot(par.asset_grid,v_star(:,1),'linewidth',2)
hold on 
plot(par.asset_grid,v_star(:,2),'linewidth',2)
legend('Low income','High income','Location','NorthWest')
xlabel('Current assets','fontsize',14) 
title('Value Function','fontsize',14)
grid on
hold off


% Plot optimal consumption
figure(2)
plot(par.asset_grid,c_pol(:,1),'linewidth',2)
hold on 
plot(par.asset_grid,c_pol(:,2),'linewidth',2)
legend('Low income','High income','Location','NorthWest')
xlabel('Current assets','fontsize',14)
ylabel('Consumption','fontsize',14) 
title('Consumption Policy Function','fontsize',14)
grid on
hold off

% Plot optimal savings with 45 degree line
figure(3)
plot(par.asset_grid,a_pol(:,1),'linewidth',2)
hold on 
plot(par.asset_grid,a_pol(:,2),'linewidth',2)
hold on 
plot(par.asset_grid,par.asset_grid,'--','linewidth',2)
legend('Low income','High income','45 line','Location','NorthWest')
xlabel('Current assets','fontsize',14)
ylabel('Next period assets','fontsize',14) 
title('Savings Policy Function','fontsize',14)
grid on
hold off

end

%% simulate a long time series and compute histogram

T = 500000;
a_sim_val = compute_asset_series(c_pol,T,par);

if do_plots==1
    figure(4)
    histogram(a_sim_val,20,'Normalization','pdf')
    title('Distribution of Asset holdings')
end

% %% How aggregate savings vary with the borrowing limit
% 
% fprintf('\n')
% fprintf('Solving Exercise 4... \n')
% 
% nr    = 25;
% r_vec = linspace(0.0,0.04,nr);
% b_vec = [1,3];
% asset_mean = nan(nr,length(b_vec));
% T      = 250000;
% max_iter = 1000;
% tol  = 1e-7; %stricter tolerance for policy iteration
% damp = 0;
% grid_max = 4;
% 
% %%% It doesn't work with b different from zero!!!
% %check if VFI works
% %re-initialize c0
% for b_c=1:length(b_vec)
%     par.b = b_vec(b_c);
%     par.asset_grid = linspace(-par.b,grid_max,na)';
%     for r_c=1:nr
%         par.R = 1+r_vec(r_c);
%         [v0,c0] = initialize(par);
%         fprintf('Doing b_c = %d and r_c = %d \n',b_c,r_c);
%         vopt   = compute_fixed_point(v0,tol,max_iter,damp,'bellman',0,par);
%         [~,cp] = bellman_operator(vopt,par);
%         asset_simul = compute_asset_series(cp,T,par);
%         asset_mean(r_c,b_c) = mean(asset_simul);
%     end
% end
% 
% %plot(par.asset_grid,cp),legend('low shock','high shock')  
% 
% figure(5)
% plot(asset_mean(:,1),r_vec,asset_mean(:,2),r_vec,'linewidth',2)
% grid on
% legend('b = 1','b = 3','location','northwest')
% xlabel('capital')
% ylabel('interest rate')
% title('How aggregate savings vary with the borrowing limit')

%% Additional functions: they are local to this script

function a_sim_val = compute_asset_series(c_optim,T,par)
%Simulates a time series of length T for assets, given optimal
%savings behavior, implied by optimal consumption

a_sim_val = zeros(T+1,1); %values
%z_sim_ind = markovchain(par.PZ,T,1,1,1); %indexes

dbg_markov = 0;
prob0V = [1,0]'; % vector of prob for init cond
random_numbers = rand(1,T);
z_sim_ind = markov_sim(1, T, prob0V, par.PZ', random_numbers, dbg_markov);

for t=1:T
    iz = z_sim_ind(t);
    a_sim_val(t+1) = par.R*a_sim_val(t)+par.z_vals(iz)-...
        myinterp1_equi(par.asset_grid,c_optim(:,iz),a_sim_val(t));
end

end %END local function

function [v0,c0] = initialize(par)

v0 = zeros(par.na,par.nz);
c0 = zeros(par.na,par.nz);
for iz=1:par.nz
    for ia=1:par.na
        c_max = par.R*par.asset_grid(ia)+par.z_vals(iz)+par.b-par.tiny;
        c0(ia,iz) = c_max;
        v0(ia,iz) = par.u(c_max)/(1-par.beta);
    end
end

end %END local function




