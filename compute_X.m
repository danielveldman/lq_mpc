function [X, duration] = compute_X(A, X0, B, U, tgrid, Mass)

% solves the forward dynamics
% Mass*\dot{x}(t) = Ax(t) + Bu(t),      x(0) = X0,
% on the time grid tgrid by the Crank-Nicholson scheme. 

tic
X = zeros(length(X0), length(tgrid));
X(:,1) = X0;
dt = diff(tgrid);
for ii = 2:length(tgrid)
    X(:,ii) = (Mass-dt(ii-1)/2*A)\(Mass*X(:,ii-1)+dt(ii-1)/2*A*X(:,ii-1) + dt(ii-1)*B*U(:,ii-1));
end
duration = toc;