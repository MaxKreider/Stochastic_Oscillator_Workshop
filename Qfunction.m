function [X, Y, Q, P0, lambda, lambda_chosen] = Qfunction(a, b, c, d, N, M, f_func, g_func, m_func, n_func)
%
% Qfunction discretizes a 2D stochastic Koopman operator with variable coefficients.
%
%   [X, Y, Q, lambda, lambda_chosen] = Qfunction(a, b, c, d, N, M, f_func, g_func, m_func, n_func)
%   builds a rectangular grid on the domain [a,b] x [c,d] with N points in the x-
%   direction and M points in the y-direction. The function then
%   discretizes the backward equation for the stochastic differential
%   equation. It also constructs the forward operator as the transpose of
%   the backward operator.
%
%   dX = m(X,Y) dt + sqrt(2*f(X,Y)) dWx(t)
%   dY = n(X,Y) dt + sqrt(2*g(X,Y)) dWy(t)
%
%   with independent Wiener process increments dWx(t) and dWy(t):
%
%   L_dagger[Q_xx] = (f_func)Q_xx + (g_func)Q_yy + (m_func)Q_x + (n_func)Q_y
%
%   at the grid points to construct finite-difference discretizations
%   (using a 9-point stencil for interior points) and adjusts for Neumann boundary conditions.
%
%   Inputs:
%       a, b   : x-domain boundaries with a < b.
%       c, d   : y-domain boundaries with c < d.
%       N, M   : number of grid points in x and y directions, respectively.
%       f_func : function handle for the coefficient f(x,y) (e.g., diffusion term in x).
%       g_func : function handle for the coefficient g(x,y) (e.g., diffusion term in y).
%       m_func : function handle for the drift term m(x,y) (affecting the x-derivative).
%       n_func : function handle for the drift term n(x,y) (affecting the y-derivative).
%
%   Outputs:
%       X, Y           : meshgrid arrays corresponding to the domain.
%       Q              : the Q-function.
%       P0             : the stationary distribution
%       lambda         : the low-lying eigenvalues of the SKO
%       lambda_chosen  : the eigenvalue corresponding to the Q-function
%
%   Example (FitzHugh-Nagumo system):
%
%       % define numerical domain:
%       a = -5;
%       b = 5;
%       c = -5;
%       d = 5;
%       N = 500;
%       M = 500;
%
%       % specify parameters for the backward equation
%       D = 0.01;
%       f_func = @(x,y) D + 0*x.*y;
%       g_func = @(x,y) D + 0*x.*y;
%       m_func = @(x,y) x - x.^3/3 - y + 1/2 + 0*x.*y;
%       n_func = @(x,y) 1/12.5*(x + 0.7 - 0.8*y) + 0*x.*y;
%
%       % generate the Q-function
%       [X, Y, Q, lambda, lambda_chosen] = Qfunction(a, b, c, d, N, M, f_func, g_func, m_func, n_func);
%
%   Author: Max Kreider
%   Date: May 8, 2025


% check the number of input arguments.
if nargin < 10
    error('Qfunction requires 10 input arguments: a, b, c, d, N, M, f_func, g_func, m_func, n_func');
end


%% create the backwards operator

% discretization
h = (b - a) / (N - 1);
k = (d - c) / (M - 1);

% create a mesh
x = linspace(a, b, N) + h/2;
y = linspace(c, d, M) + k/2;
[X, Y] = meshgrid(x, y);

% reshape for function evaluation
X1 = reshape(X', N*M, 1);
Y1 = reshape(Y', N*M, 1);

% evalute variable coefficients
V1 = f_func(X1, Y1);
V2 = g_func(X1, Y1);
V3 = m_func(X1, Y1);
V4 = n_func(X1, Y1);

% coefficients for the boundary
coeff_A = 2*(V1 + V2*(h/k)^2);
coeff_B = -V1 + 1/2*h*V3;
coeff_C = -V1 - 1/2*h*V3;
coeff_D = -(h/k)^2*V2 - h^2/(2*k)*V4;
coeff_E = -(h/k)^2*V2 + h^2/(2*k)*V4;

% coefficients for the inner points
coeff_a = 1/12*V1 - h/12*V3;
coeff_b = -16/12*V1 + 8/12*h*V3;
coeff_c = -16/12*V1 - 8/12*h*V3;
coeff_d = 1/12*V1 + h/12*V3;
coeff_i = 1/12*(h/k)^2*V2 - 1/12*(h^2/k)*V4;
coeff_h = -16/12*(h/k)^2*V2 + 8/12*(h^2/k)*V4;
coeff_f = -16/12*(h/k)^2*V2 - 8/12*(h^2/k)*V4;
coeff_e = 1/12*(h/k)^2*V2 + 1/12*(h^2/k)*V4;
coeff_j = 30/12*(V1 + (h/k)^2*V2);

% adjust boundary coefficients due to Neumann condition
coeff_D(1:N) = coeff_D(1:N) + coeff_E(1:N);
coeff_E(end-(N-1):end) = coeff_E(end-(N-1):end) + coeff_D(end-(N-1):end);
coeff_B(N:N:end) = coeff_B(N:N:end) + coeff_C(N:N:end);
coeff_C(1:N:end) = coeff_C(1:N:end) + coeff_B(1:N:end);

% construct the outer matrix (to account for boundary conditions)
e = ones(N*M, 1);
e(2*N+1:(M-2)*N) = 0;
e(2*N+1:N:(M-2)*N) = 1;
e(3*N:N:(M-2)*N) = 1;
e(2*N+2:N:(M-2)*N) = 1;
e(3*N-1:N:(M-2)*N) = 1;
e2 = e;
e3 = e;
e2(N:N:end) = 0;
e3(N+1:N:end) = 0;

% discretized version of backward operator -- the part from the boundary:
A_out = spdiags([coeff_D.*e, coeff_C.*e2, coeff_A.*e, coeff_B.*e3, coeff_E.*e], ...
    [-N, -1, 0, 1, N], N*M, N*M)';

% construct the inner matrix
e = zeros(N*M, 1);
e(2*N+3:(M-2)*N-2) = 1;
e(2*N+1:N:(M-2)*N) = 0;
e(3*N:N:(M-2)*N) = 0;
e(2*N+2:N:(M-2)*N) = 0;
e(3*N-1:N:(M-2)*N) = 0;
e2 = e;
e3 = e;
e4 = e;
e5 = e;
e2(N:N:end) = 0;
e3(N+1:N:end) = 0;
e4(N-1:N:end) = 0;
e5(N+2:N:end) = 0;

% discretized version of backward operator -- the part from the interior:
A_in = spdiags([coeff_e.*e, coeff_f.*e, coeff_d.*e4, coeff_c.*e2, coeff_j.*e, ...
    coeff_b.*e3, coeff_a.*e5, coeff_h.*e, coeff_i.*e], ...
    [-2*N, -N, -2, -1, 0, 1, 2, N, 2*N], N*M, N*M)';

% form the backward operator
A_b = -1/h^2*(A_in + A_out);

% form forward operator
A_f = A_b';


%% diagonalize the operators

% diagonalize
[VV_b, lambda] = eigs(A_b,15,1e-4);
[VV_f,lambda_f] = eigs(A_f,15,1e-4);

% normalize
Q_b = VV_b./max(VV_b, [], 1);
Q_f = VV_f./max(VV_f, [], 1);

% sort the Q eigenvalues (only positive imaginary part)
lambda = diag(lambda);
imag_index = find(imag(lambda)>0);
lambda_temp = lambda(imag_index);
Q_b = Q_b(:,imag_index);

% sort the Q eigenvalues (keep only the one with least negative real part)
real_index = find(max(real(lambda_temp)));
lambda_chosen = lambda_temp(real_index);
Q = Q_b(:,real_index);

% reshape Q
Q = reshape(Q,[N,M]).';

% sort the forward eigenvalues
lambda_f = diag(lambda_f);
[~, stat_index] = min(abs(lambda_f));
P0 = Q_f(:,stat_index);

% reshape P0
P0 = reshape(P0,[N,M]).';
P0 = P0/sum(sum(P0)*h*k);

% normalize Q
I = trapz(y,trapz(x,abs(Q).^2.*P0,2));
Q = Q/sqrt(I);

end
