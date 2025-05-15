function[t,u] = traj_StuartLandau(D)
%
% function [t,u]=traj_StuartLandau(D);
%
% This function numerically integrates a (noisy) FitzHugh-Nagumo oscillator:
%
%   dx = [-4x(x^2 + y^2 - 1) + 2y]dt + sqrt(2D)dW_1(t)
%   dy = [-4y(x^2 + y^2 - 1) - 2x]dt + sqrt(2D)dW_2(t)
%
% Dependencies:
%   - TimeSeries.m (for SDE simulation)
%
% Usage:
%   Input D = 0 for deterministic system (default value if no input given)
%   Input D = 10^-4 for small noise
%   Input D = 10^-3 for medium noise
%   Input D = 10^-2 for larger noise
%
% Figures:
%   - Figure 1 displays time-series data in the original coordinates
%
% Example:
%   D = 0
%   traj_StuartLandau(D);
%
% Author: Max Kreider
% Date: May 8, 2025


%% check user input

%check user input
if nargin == 0
    D = 0;
elseif nargin == 1 && D<0
    fprintf('\n\n Please enter D>=0 \n\n')
end


%% generate time series

%display progress update
fprintf('\n\nGenerating time-series data in original coordinates... \n\n')

%define numerical domain
a = -3;
b = 3;
c = -3;
d = 3;

%drift and diffusion terms
f = @(t,y)[-4*y(1)*(y(1)^2+y(2)^2-1)+2*y(2); -4*y(2)*(y(1)^2+y(2)^2-1)-2*y(1)];
g = @(t,y)[sqrt(2*D); sqrt(2*D)];

%simulation parameters
tmax = 50;
dt = 1/256;
y0 = [1; 0.1];

%run the simulation
[t, u] = TimeSeries(f, g, tmax, dt, y0, 'BC', a, b, c, d);

%for computing nullclines
m_func = @(x,y) -4*x.*(x.^2+y.^2-1)+2*y + 0*x.*y;
n_func = @(x,y) -4*y.*(x.^2+y.^2-1)-2*x + 0*x.*y;

%domain
N = 400;
M = 400;
h = (b - a) / (N - 1);
k = (d - c) / (M - 1);
x = linspace(a, b, N) + h/2;
y = linspace(c, d, M) + k/2;
[X, Y] = meshgrid(x, y);


%% visualize (if needed)

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
axis equal tight
grid on
set(gca,'FontSize',15)
box on

end
