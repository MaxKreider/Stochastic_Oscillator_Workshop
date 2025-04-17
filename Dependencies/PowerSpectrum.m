function[power_x,power_y,power_Q,power_exact_Q] = PowerSpectrum(f, g, Tmax, dt, Num, freq, M, y0, X, Y, Q, lambda)


%% compute the power spectra from time series

% initialize power spectra in original coordinates, and Q-function coordinates
power_x = 0;
power_y = 0;
power_Q = 0;

% compute 
for i=1:M
    [~, u] = TimeSeries(f, g, Tmax, dt, y0);
    power_x = power_x + abs(fftshift(fft(u(1,:)))).^2;
    power_y = power_y + abs(fftshift(fft(u(2,:)))).^2;
    power_Q = power_Q + abs(fftshift(fft(interp2(X, Y, Q, u(1,:), u(2,:), 'linear')))).^2;
end

% normalize
power_x = power_x/M/Num*dt;
power_y = power_y/M/Num*dt;
power_Q = power_Q/M/Num*dt;


%% compute exact power spectrum in Q-function coordinates

%exact solution
mu = real(lambda);
omega = imag(lambda);
power_exact_Q = 2*abs(mu)./(mu^2+(freq-omega).^2);

end