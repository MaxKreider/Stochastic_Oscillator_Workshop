function 
% This script demonstrates the following steps for a (noisy) FitzHugh-Nagumo oscillator:
%
%   dx = [x - x^3/3 - y + 1/2]dt + sqrt(2D)dW_1(t)
%   dy = [(x + 0.7 - 0.8*y)/12.5]dt + sqrt(2D)dW_2(t) 
%
%   1. Simulating a stochastic differential equation (SDE) with the function TimeSeries.m
%
% Dependencies:
%   - TimeSeries.m (for SDE simulation)
%
% Usage:
%   Simply run the script to execute the full workflow. Adjust parameters for different models as needed.
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
a = -2.5; 
b = 2.5; 
c = -2.5; 
d = 2.5;

%parameters
%D = 0.1215;     %large noise
%D = 0.05;     %medium noise
%D = 0.0005;     %small noise
D = 0;          %no noise


%drift and diffusion terms
f = @(t,y)[y(1)-y(1).^3/3 - y(2) + 0.5; (y(1) + 0.7  - 0.8*y(2))/12.5];
g = @(t,y)[sqrt(2*D); sqrt(2*D)];

%simulation parameters
tmax = 500;
dt = 1/256;
y0 = [1; 0.1];

%run the simulation
[t, u] = TimeSeries(f, g, tmax, dt, y0, 'BC', a, b, c, d);

%for computing nullclines
m_func = @(x,y) x - x.^3/3 - y + 1/2 + 0*x.*y;
n_func = @(x,y) 1/12.5*(x + 0.7 - 0.8*y) + 0*x.*y;

%domain
N = 400;
M = 400;
h = (b - a) / (N - 1);
k = (d - c) / (M - 1);
x = linspace(a, b, N) + h/2;
y = linspace(c, d, M) + k/2;
[X, Y] = meshgrid(x, y);


%% visualize (if needed)

%display progress update
fprintf('Generating plots, if requested by user input ... \n\n')
fprintf('................................................ \n\n')

%time series in original coordinates
reply = input('Display time series in original coordinates? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')

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