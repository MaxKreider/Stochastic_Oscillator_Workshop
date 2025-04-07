%
% This script demonstrates the following steps for a noisy Morris-Lecar neuron:
%
%   dx = [1/C*(I-gL*(x-vL) - gK*y*(x-vK) - gCA*m_inf(x)*(x-vCA))]dt + sqrt(2Dv)dW_1(t)
%   dy = [alpha(x)*(1-y) - beta(x)*y]dt + [sqrt(2*Dn^2/2*(alpha(x)*(1-y)+beta(x)*y))]dW_2(t) 
%
%   1. Simulating a stochastic differential equation (SDE) with the function TimeSeries.m
%   2. Constructing the Q-function by discretizing the SKO on a specified domain with the function Qfunction.m
%   3. Interpolating the SDE time series onto the Q-function coordinate system
%   4. Visualizing the time series in both the original and Q-function coordinates, as well as the SKO spectrum 
%      and the Q-function
%
% Dependencies:
%   - TimeSeries.m (for SDE simulation)
%   - Qfunction.m (for constructing the Q-function)
%
% Usage:
%   Simply run the script to execute the full workflow. Adjust parameters for different models as needed.
%
% Figures:
%   - Figure 1 displays time-series data (top: x(t), bottom y(t)) in the original coordinates
%   - Figure 2 displays time-series data (top: Re(Q(t)), bottom Im(Q(t)) in the Q-function coordinates
%   - Figure 3 displays the low-lying eigenvalue spectrum of the SKO. 
%       The Q-function eigenvalue is highlighted in pink
%   - Figures 4 and 5 display the real and imaginary parts of the Q-function
%   - Figure 6 displays the asymptotic stochastic phase
%
% Author: Max Kreider
% Date: April 7, 2025


%% generate time series

% parameter values
global I Dn vK vL vCA gK gL gCA vA vB vC vD C phi Dv

I = 180;
Dn = 1e-1;
Dv = 1;
vK = -84;
vL = -60;
vCA = 120;
gK = 8;
gL = 2;
gCA = 4.4;
vA = -1.2;
vB = 18;
vC = 2;
vD = 30;
C = 20;
phi = 0.04;

% drift and diffusion terms
f = @(t,y)[1/C*(I-gL*(y(1)-vL)-gK*y(2).*(y(1)-vK)-gCA.*m_inf(y(1),vA,vB).*(y(1)-vCA)); (alpha(y(1),phi,vC,vD).*(1-y(2))-beta(y(1),phi,vC,vD).*y(2))];
g = @(t,y)[sqrt(2*Dv); sqrt(2*Dn^2/2*(alpha(y(1),phi,vC,vD).*(1-y(2))+beta(y(1),phi,vC,vD).*y(2)))];

% simulation parameters
tmax = 1000;
dt = 1/256;
y0 = [1; 0.5];

% run the simulation
[t, u] = TimeSeries(f, g, tmax, dt, y0);


%% construct the Q-function

% define numerical domain
a = -70;
b = 60;
c = 0;
d = 1;
N = 1000;
M = 400;

% specify parameters for the backward equation
f_func = @(x,y) Dv + 0*x.*y;
g_func = @(x,y) Dn^2/2*(alpha(x,phi,vC,vD).*(1-y)+beta(x,phi,vC,vD).*y) + 0*x.*y;
m_func = @(x,y) 1/C*(I-gL*(x-vL)-gK*y.*(x-vK)-gCA.*m_inf(x,vA,vB).*(x-vCA)) + 0*x.*y;
n_func = @(x,y) (alpha(x,phi,vC,vD).*(1-y)-beta(x,phi,vC,vD).*y) + 0*x.*y;

% generate the Q-function
[X, Y, Q, lambda, lambda_chosen] = Qfunction(a, b, c, d, N, M, f_func, g_func, m_func, n_func);


%% interpolate time series in Q-function coordinates

%time series in Q-function coordinates
Q_series = interp2(X, Y, Q, u(1,:), u(2,:), 'linear');


%% visualize (if needed)

% time series in original coordinates
figure(1)
subplot(2,1,1)
plot(t,u(1,:),'k','LineWidth',2)
xlabel('time t')
ylabel('x(t)')
title('Time series, original coordinates')
set(gca,'fontsize',15)
subplot(2,1,2)
plot(t,u(2,:),'k','LineWidth',2)
xlabel('time t')
ylabel('y(t)')
set(gca,'fontsize',15)

% time series in Q-function coordinates
figure(2)
subplot(2,1,1)
plot(t,real(Q_series),'k','LineWidth',2)
xlabel('time t')
ylabel('Re(Q(t))')
title('Time series, Q-function coordinates')
set(gca,'fontsize',15)
subplot(2,1,2)
plot(t,imag(Q_series),'k','LineWidth',2)
xlabel('time t')
ylabel('Im(Q(t))')
set(gca,'fontsize',15)

% low-lying SKO eigenvalues (and highlight Q-function eigenvalue in pink)
figure(3)
hold on
plot(real(lambda),imag(lambda),'k.','MarkerSize',40)
plot(real(lambda_chosen), imag(lambda_chosen), 'm.', 'MarkerSize', 30)
plot(real(lambda_chosen), imag(-lambda_chosen), 'm.', 'MarkerSize', 30)
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title('Spectrum')
box on
set(gca,'FontSize',15)

% Re(Q)
figure(4)
contourf(X,Y,real(Q),500,'LineColor','none')
colormap jet
colorbar
xlabel('x')
ylabel('y')
title('Re(Q(x,y))')
box on
axis square
set(gca,'FontSize',15)

% Im(Q)
figure(5)
contourf(X,Y,imag(Q),500,'LineColor','none')
colormap jet
colorbar
xlabel('x')
ylabel('y')
title('Im(Q(x,y))')
box on
axis square
set(gca,'FontSize',15)

% Arg(Q)
figure(6)
contourf(X,Y,angle(Q)+pi,500,'LineColor','none')
colormap jet
colorbar
xlabel('x')
ylabel('y')
title('Arg(Q(x,y))')
box on
axis square
set(gca,'FontSize',15)


%% functions

%alpha
function[u]=alpha(V,phi,vC,vD)
xi=(V-vC)/vD;
u=phi*cosh(xi/2)./(1+exp(-2*xi));
end

%beta
function[u]=beta(V,phi,vC,vD)
xi=(V-vC)/vD;
u=phi*cosh(xi/2)./(1+exp(2*xi));
end

%m
function[u]=m_inf(V,vA,vB)
xi=(V-vA)/vB;
u=1/2*(1+tanh(xi));
end

 