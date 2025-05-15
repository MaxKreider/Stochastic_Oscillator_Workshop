function [t,u]=traj_MorrisLecar(Dv_input,Dn_input)

%function [t,u]=traj_MorrisLecar(Dv,Dn)
%
% This function numerically integrates a (noisy) Morris-Lecar neuron:
%
%   dx = [1/C*(I-gL*(x-vL) - gK*y*(x-vK) - gCA*m_inf(x)*(x-vCA))]dt + sqrt(2Dv)dW_1(t)
%   dy = [alpha(x)*(1-y) - beta(x)*y]dt + [sqrt(2*Dn^2/2*(alpha(x)*(1-y)+beta(x)*y))]dW_2(t)
%
% Dependencies:
%   - TimeSeries.m (for SDE simulation)
%
% Usage:
%       % Choose a noise setting
%       % Dn=0; Dv=0; No noise (default)
%       % Dn = 5*1e-2; Dv = .5; (small noise in both n-gate and voltage components)
%       % Dn = 8*1e-1; Dv = 2; (large noise in both n-gate and voltage components)
%       [t,u]=traj_MorrisLecar(Dv,Dn)
%
% Figures:
%   - Figure 1 displays time-series data in the original coordinates
%
% Author: Max Kreider
% Date: May 8, 2025

%% generate time series

%display progress update
fprintf('\n\nGenerating time-series data in original coordinates... \n\n')

%define numerical domain
a = -90;
b = 90;
c = 0;
d = 1;

%parameter values
global I Dn vK vL vCA gK gL gCA vA vB vC vD C phi Dv

if nargin > 0, Dv=Dv_input; else Dv=0; end % default to zero noise case
if nargin > 1, Dn=Dn_input; else Dn=0; end % default to zero noise case

I = 180;
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

%drift and diffusion terms
f = @(t,y)[1/C*(I-gL*(y(1)-vL)-gK*y(2).*(y(1)-vK)-gCA.*m_inf(y(1),vA,vB).*(y(1)-vCA)); (alpha(y(1),phi,vC,vD).*(1-y(2))-beta(y(1),phi,vC,vD).*y(2))];
g = @(t,y)[sqrt(2*Dv); sqrt(2*Dn^2/2*(alpha(y(1),phi,vC,vD).*(1-y(2))+beta(y(1),phi,vC,vD).*y(2)))];

%simulation parameters
tmax = 1000;
dt = 1/256;
y0 = [1; 0.5];

%run the simulation
[t, u] = TimeSeries(f, g, tmax, dt, y0, 'BC', a, b, c, d);

%for computing nullclines
m_func = @(x,y) 1/C*(I-gL*(x-vL)-gK*y.*(x-vK)-gCA.*m_inf(x,vA,vB).*(x-vCA)) + 0*x.*y;
n_func = @(x,y) (alpha(x,phi,vC,vD).*(1-y)-beta(x,phi,vC,vD).*y) + 0*x.*y;

%domain
N = 500;
M = 500;
h = (b - a) / (N - 1);
k = (d - c) / (M - 1);
x = linspace(a, b, N) + h/2;
y = linspace(c, d, M) + k/2;
[X, Y] = meshgrid(x, y);


%% visualize (if needed)

%time series in original coordinates
figure(1)
set(gcf,'position',[66.60000000000001,163.4,899.2,420])

%Left column: two stacked subplots for x(t) and y(t)
subplot(2,2,1)
plot(t, u(1,:), 'k', 'LineWidth', 2)
ylabel('x(t)')
title('x(t) vs t')
xlim([0 tmax])
set(gca,'FontSize',15)

subplot(2,2,3)
plot(t, u(2,:), 'k', 'LineWidth', 2)
xlabel('time t')
ylabel('y(t)')
title('y(t) vs t')
xlim([0 tmax])
set(gca,'FontSize',15)

%Right column: phase‐plane trajectory spanning both rows
subplot(2,2,[2 4])
hold on
plot(u(1,:), u(2,:), 'k', 'LineWidth', 2) %trajectories
contour(X, Y, m_func(X,Y), [0 0], 'g', 'LineWidth', 2);  % x-nullcline in green
contour(X, Y, n_func(X,Y), [0 0], 'm', 'LineWidth', 2);  % y-nullcline in pink
xlabel('x')
ylabel('y')
title('Phase‐Plane and Nullclines')
%axis equal tight
grid on
set(gca,'FontSize',15)
box on


end

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