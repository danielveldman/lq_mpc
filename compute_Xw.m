function [X, duration] = compute_Xw(A, X0, B, U, tgrid, Mass, w)

% solves the forward dynamics
% Mass*\dot{x}(t) = Ax(t) + Bu(t) + w(t),      x(0) = X0,
% on the time grid tgrid by the Crank-Nicholson scheme. 

tic
X = zeros(length(X0), length(tgrid));
X(:,1) = X0;
dt = diff(tgrid);
for ii = 2:length(tgrid)
    X(:,ii) = (Mass-dt(ii-1)/2*A)\(Mass*X(:,ii-1)+dt(ii-1)/2*A*X(:,ii-1) + dt(ii-1)*B*U(:,ii-1) + dt(ii-1)*(w(tgrid(ii-1))+w(tgrid(ii))/2));
end
duration = toc;