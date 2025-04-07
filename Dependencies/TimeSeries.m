function [t, y] = TimeSeries(f, g, tmax, dt, y0)
% TimeSeries simulates an SDE using the Euler-Maruyama method.
%
%   [t, y] = TimeSeries(f, g, tmax, dt, y0) integrates the stochastic
%   differential equation
%
%   dy = f(x,t)dt + g(x,t)dW(t)
%
%   defined by the drift function f and the diffusion
%   function g over the interval [0, tmax] with time step dt and initial
%   condition y0.
%
%   Inputs:
%       f    - function handle for the drift term, f(t,y)
%       g    - function handle for the diffusion term, g(t,y)
%       tmax - maximum simulation time
%       dt   - time step size
%       y0   - column vector of initial conditions
%
%   Outputs:
%       t    - time vector
%       y    - solution matrix; each column corresponds to the solution
%              at the corresponding time in t.
%
%   Example (FitzHugh-Nagumo system):
%
%       % Define drift and diffusion functions
%       f = @(t,y)[y(1) - y(1).^3/3 - y(2) + 0.5; (y(1) + 0.7 - 0.8*y(2))/12.5];
%       D = 0.01; % diffusion coefficient used inside g if needed
%       g = @(t,y)[sqrt(2*D); sqrt(2*D)];
%
%       % Set simulation parameters
%       tmax = 1000;
%       dt = 1/256;
%       y0 = [1; 0];
%
%       % Run simulation
%       [t, y] = TimeSeries(f, g, tmax, dt, y0);
%
%   Author: Max Kreider
%   Date: April 7, 2025

% check the number of input arguments.
if nargin < 5
    error('TimeSeries requires 5 input arguments: f, g, tmax, dt, y0');
end

% create time vector
t = 0:dt:tmax;

% extract problem dimensions
numSteps = length(t);
numVars = length(y0);

% preallocate solution matrix
y = zeros(numVars, numSteps);
y(:,1) = y0;

% Euler-Maruyama integration
for k = 1:numSteps-1

    % evaluate drift and diffusion at current time and state
    drift = f(t(k), y(:,k));
    diffusion = g(t(k), y(:,k));

    % update the state vector
    y(:,k+1) = y(:,k) + dt*drift + sqrt(dt)*diffusion.*randn(numVars, 1);
end

end