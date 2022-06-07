function [phi, duration] = compute_phi(A, Q, X, xd, tgrid, Mass)

% solves the backward dynamics
% Mass.'*\dot{phi}(t) = A.'phi(t) + Q*(x(t) - xd(t)),      phi(T) = 0,
% with a scheme that leads to discretely consistent gradients when the
% Crank-Nicholson scheme is used for the forward dynamics.

% see:
% Apel and Flaig, Crank--Nicolson Schemes for Optimal Control Problems
% with Evolution Equations, SIAM journal on Numerical Analysis, 50(3), 2012. 

tic
dt = diff(tgrid); ndt = length(dt);
phi = zeros(size(X,1), length(dt)); N = size(X,1);
dx = X - xd(tgrid.'); I = speye(N);

phi(:,end) = (I-dt(end)/2*A.')\((dt(end)*Q*dx(:,end)/2));
for ii = ndt-1:-1:1
    phi(:,ii) = (Mass-dt(ii)/2*A.')\(Mass*phi(:,ii+1)+dt(ii)/2*A.'*phi(:,ii+1) + dt(ii)*Q*dx(:,ii+1));
end
duration = toc;