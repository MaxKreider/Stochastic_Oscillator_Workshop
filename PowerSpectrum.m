function[power_x,power_y,power_Q,power_exact_Q] = PowerSpectrum(f, g, Tmax, dt, Num, freq, M, y0, X, Y, Q, lambda, flag, a, b, c, d)
%
% PowerSpectrum computes the power spectrum of a stochastic oscillator.
%
%   [power_x, power_y, power_Q, power_exact_Q] = PowerSpectrum(f, g, Tmax, dt, Num, freq, M, y0, X, Y, Q, lambda, flag, a, b, c, d)
%   simulates a stochastic differential equation multiple times and computes
%   the power spectra of the time series in the original coordinates (x and y)
%   and in transformed coordinates using the Q-function, Q(x,y). The exact theoretical power spectrum in 
%   Q-function coordinates is also returned.
%
%   Inputs:
%       f           - drift function handle, f(t,y)
%       g           - diffusion function handle, g(t,y)
%       Tmax        - maximum simulation time
%       dt          - time step size
%       Num         - number of time steps in each simulation
%       freq        - frequency vector for evaluating the exact power spectrum
%       M           - number of Monte Carlo trials
%       y0          - initial condition for the stochastic simulation
%       X, Y        - meshgrid arrays corresponding to the domain of Q
%       Q           - 2D array defining Q(X,Y), used to transform coordinates
%       lambda      - complex eigenvalue associated with Q, used in the exact spectrum
%       flag        - optional string flag, set to 'BC' to enable 
%                     reflecting boundary conditions (default is no boundaries)
%       a,b,c,d     - specifies domain [a,b] x [c,d] on which to implement
%                     reflecting boundary conditions
%
%   Outputs:
%       power_x         - empirical power spectrum for x(t)
%       power_y         - empirical power spectrum for y(t)
%       power_Q         - empirical power spectrum for Q(x(t),y(t))
%       power_exact_Q   - theoretical power spectrum for Q
%
%   Dependencies:
%       - TimeSeries.m (for SDE simulation)
%       - Qfunction.m (for constructing the Q-function)
%
%   Example (FitzHugh-Nagumo system):
%
%       % parameters
%       D = 0.0005;
%
%       % drift and diffusion terms
%       f = @(t,y)[y(1)-y(1).^3/3 - y(2) + 0.5; (y(1) + 0.7  - 0.8*y(2))/12.5];
%       g = @(t,y)[sqrt(2*D); sqrt(2*D)];
%
%       % initial condition
%       y0 = [1; 0.1];
%
%       % define numerical domain
%       a = -3; 
%       b = 3; 
%       c = -3; 
%       d = 3;
%       N = 400; 
%       M = 400;
%
%       % specify parameters for the backward equation
%       f_func = @(x,y) D + 0*x.*y;
%       g_func = @(x,y) D + 0*x.*y;
%       m_func = @(x,y) x - x.^3/3 - y + 1/2 + 0*x.*y;
%       n_func = @(x,y) 1/12.5*(x + 0.7 - 0.8*y) + 0*x.*y;
%
%       % generate the Q-function
%       [X, Y, Q, P0, lambda, lambda_chosen] = Qfunction(a, b, c, d, N, M, f_func, g_func, m_func, n_func);
%
%       % time
%       Delta = 1/50;
%       Num = 2^18;
%       pst = 0:Delta:(Num-1)*Delta;
%
%       % frequency vector
%       step = (-Num/2:Num/2-1);
%       freq = 1/(Num*Delta)*step*2*pi;
%
%       % number of trials
%       M = 100;
%
%       % compute power spectra
%       [power_x,power_y,power_Q,power_exact_Q] = PowerSpectrum(f, g, pst(end), Delta, Num, freq, M, y0*rand, X, Y, Q, lambda_chosen);
%
%   Author: Max Kreider
%   Date: May 8, 2025


% check the number of input arguments
if nargin < 12
    error('PowerSpectrum.m requires at least 12 input arguments: f, g, Tmax, dt, Num, freq, M, y0, X, Y, Q, lambda');
elseif nargin == 12
    flag = 0;  % no boundaries given
elseif nargin == 17
    if ischar(flag) && strcmp(flag, 'BC')
        flag = 1;  % boundary values provided
    else
        error('PowerSpectrum.m only supports flag = BC');
    end
else
    error('PowerSpectrum.m requires either 12 or 17 input arguments.');
end


%% compute the power spectra from time series

% initialize power spectra in original coordinates, and Q-function coordinates
power_x = 0;
power_y = 0;
power_Q = 0;

% compute
if flag == 1
    for i=1:M
        [~, u] = TimeSeries(f, g, Tmax, dt, y0, 'BC', a, b, c, d);
        power_x = power_x + abs(fftshift(fft(u(1,:)))).^2;
        power_y = power_y + abs(fftshift(fft(u(2,:)))).^2;
        power_Q = power_Q + abs(fftshift(fft(interp2(X, Y, Q, u(1,:), u(2,:), 'spline')))).^2;
    end
else
    for i=1:M
        [~, u] = TimeSeries(f, g, Tmax, dt, y0);
        power_x = power_x + abs(fftshift(fft(u(1,:)))).^2;
        power_y = power_y + abs(fftshift(fft(u(2,:)))).^2;
        power_Q = power_Q + abs(fftshift(fft(interp2(X, Y, Q, u(1,:), u(2,:), 'spline')))).^2;
    end
end

% normalize
power_x = power_x/M/Num*dt;
power_y = power_y/M/Num*dt;
power_Q = power_Q/M/Num*dt;


%% compute exact power spectrum in Q-function coordinates

% exact solution
mu = real(lambda);
omega = imag(lambda);
power_exact_Q = 2*abs(mu)./(mu^2+(freq-omega).^2);

end
