%
% This script demonstrates the following steps for a noisy Stuart-Landau oscillator:
%
%   dx = [-4x(x^2 + y^2 - 1) + 2y]dt + sqrt(2D)dW_1(t)
%   dy = [-4y(x^2 + y^2 - 1) - 2x]dt + sqrt(2D)dW_2(t)
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
% Date: April 17, 2025


%% generate time series

%display progress update
fprintf('\n\nGenerating time-series data in original coordinates... \n\n')

% parameters
D = 0.09473;

% drift and diffusion terms
f = @(t,y)[-4*y(1)*(y(1)^2+y(2)^2-1)+2*y(2); -4*y(2)*(y(1)^2+y(2)^2-1)-2*y(1)];
g = @(t,y)[sqrt(2*D); sqrt(2*D)];

% simulation parameters
tmax = 50;
dt = 1/256;
y0 = [1; 0];

% run the simulation
[t, u] = TimeSeries(f, g, tmax, dt, y0);


%% construct the Q-function

%display progress update
fprintf('Generating the Q-function and the low-lying SKO eigenspectrum... \n\n')

% define numerical domain
a = -3;
b = 3;
c = -3;
d = 3;
N = 400;
M = 400;

% specify parameters for the backward equation
f_func = @(x,y) D + 0*x.*y;
g_func = @(x,y) D + 0*x.*y;
m_func = @(x,y) -4*x.*(x.^2+y.^2-1)+2*y + 0*x.*y;
n_func = @(x,y) -4*y.*(x.^2+y.^2-1)-2*x + 0*x.*y;

% generate the Q-function
[X, Y, Q, P0, lambda, lambda_chosen] = Qfunction(a, b, c, d, N, M, f_func, g_func, m_func, n_func);


%% interpolate time series in Q-function coordinates

%display progress update
fprintf('Generating time-series in Q-function coordinates... \n\n')

% time series in Q-function coordinates
Q_series = interp2(X, Y, Q, u(1,:), u(2,:), 'linear');


%% power spectra

%display progress update
fprintf('Generating power spectra... \n\n')

% time
Delta = 1/100;
Num = 2^17;
pst = 0:Delta:(Num-1)*Delta;

% frequency vector
step = (-Num/2:Num/2-1);
freq = 1/(Num*Delta)*step*2*pi;

% number of trials
M = 100;

% compute power spectra
[power_x,power_y,power_Q,power_exact_Q] = PowerSpectrum(f, g, pst(end), Delta, Num, freq, M, y0*rand, X, Y, Q, lambda_chosen);


%% visualize (if needed)

%display progress update
fprintf('Generating plots, if requested by user input ... \n\n')
fprintf('................................................ \n\n')

% time series in original coordinates
reply = input('Display time series in original coordinates? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')

    figure(1)
    set(gcf,'position',[66.60000000000001,163.4,899.2,420])

    % Left column: two stacked subplots for x(t) and y(t)
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

    % Right column: phase‐plane trajectory spanning both rows
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

% time series in Q-function coordinates
fprintf('\n\n')
reply = input('Display time series in Q-function coordinates? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')
    
    figure(2)
    set(gcf,'position',[66.60000000000001,163.4,899.2,420])

    % Left column: two stacked subplots for x(t) and y(t)
    subplot(2,2,1)
    plot(t, real(Q_series), 'k', 'LineWidth', 2)
    ylabel('Re(Q(t))')
    title('Time series in Q-function coordinates')
    xlim([0 tmax])
    set(gca,'FontSize',15)

    subplot(2,2,3)
    plot(t, imag(Q_series), 'k', 'LineWidth', 2)
    xlabel('time t')
    ylabel('Im(Q(t))')
    xlim([0 tmax])
    set(gca,'FontSize',15)

    % Right column: phase‐plane trajectory spanning both rows
    subplot(2,2,[2 4])
    hold on
    plot(real(Q_series), imag(Q_series), 'k', 'LineWidth', 2) %trajectories
    %contour(X, Y, m_func(X,Y), [0 0], 'g', 'LineWidth', 2);  % x-nullcline in green
    %contour(X, Y, n_func(X,Y), [0 0], 'm', 'LineWidth', 2);  % y-nullcline in pink
    xlabel('Re(Q(t))')
    ylabel('Im(Q(t))')
    title('Real and imaginary time series')
    axis equal tight
    grid on
    set(gca,'FontSize',15)
    box on
end

% plot power spectra
fprintf('\n\n')
reply = input('Display power spectra? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')
    figure(3)
    hold on
    plot(freq,power_x,'-','color',[0.6 0.6 0.6 .7],'linewidth',5)
    plot(freq,power_y,'-','color',[0.6 0.6 0.6 .7],'linewidth',5)
    plot(freq,power_Q,'-','color',[0.4940 0.1840 0.5560 1],'linewidth',10)
    plot(freq,power_exact_Q,'-','color','m','linewidth',2)
    xlim([imag(lambda_chosen)-imag(lambda_chosen)*.8 imag(lambda_chosen)+imag(lambda_chosen)*.8])
    ylim([0 3*max(power_exact_Q)])
    xlabel('frequency \nu')
    ylabel('S_1(\nu)')
    title('Power spectra')
    box on
    axis square
    set(gca,'fontsize',15)
    legend('original coordinates','','Q-function coordinates','analytic expression')
end

% low-lying SKO eigenvalues (and highlight Q-function eigenvalue in pink)
fprintf('\n\n')
reply = input('Display eigenvalues? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')
    figure(4)
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
end

% Re(Q)
fprintf('\n\n')
reply = input('Display the Q-function? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')
    figure(5)
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
    figure(6)
    contourf(X,Y,imag(Q),500,'LineColor','none')
    colormap jet
    colorbar
    xlabel('x')
    ylabel('y')
    title('Im(Q(x,y))')
    box on
    axis square
    set(gca,'FontSize',15)
end

% Arg(Q)
fprintf('\n\n')
reply = input('Display the stochastic asymptotic phase? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')
    figure(7)
    contourf(X,Y,angle(Q)+pi,500,'LineColor','none')
    colormap jet
    colorbar
    xlabel('x')
    ylabel('y')
    title('Arg(Q(x,y))')
    box on
    axis square
    set(gca,'FontSize',15)
end

% Stationary distribution
fprintf('\n\n')
reply = input('Display the stationary distribution and 10 isochrons? (y = yes, any other key = no): ','s');
if strcmpi(reply,'y')
    fig = figure(8);
    ax1 = axes(fig);
    ax2 = copyobj(ax1,fig);
    
    contourf(ax1,X,Y,P0,500,'LineColor','none')
    contour(ax2,X,Y,angle(Q)+pi,10,'w','LineWidth',2.2)

    colormap(ax1,'jet')
    colormap(ax2,'gray')

    ax2.UserData = linkprop([ax1,ax2],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
    ax2.Visible = 'off';

    colorbar
    xlabel('x')
    ylabel('y')
    title('Isochrons and Stationary Distribution')
    box on
    set(gca,'FontSize',15)
    axis(ax1,'square')
    axis(ax2,'square')
end
