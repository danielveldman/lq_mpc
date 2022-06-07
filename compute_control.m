function [u0, J0, duration] = compute_control(A, X0, B, u0, Q, R, xd, tgrid, Mass)

% a code that computes the optimal control u0 that minimizes the functional
% \int (x(t) - xd(t))^T Q (x(t)-xd(t)) + u(t)^T R u(t) dt
% subject to the dynamics
% Mass*\dot{x}(t) = Ax(t) + Bu(t),      x(0) = X0,
% on the time grid tgrid. 
% The input u0 is the initial guess for the control. 
% The control is computed by a steepest descent algorithm in which the step
% size is chosen optimally. 

% This function calls the external functions:
% 1) compute_X (to solve the forward dynamics)
% 2) compute_phi (to solve the backward dynamics/adjoint state)

dt = diff(tgrid);
tic

max_iters = 1000; tol = 1e-3;
for ii = 1:max_iters
    xsol0 = compute_X(A, X0, B, u0, tgrid, Mass);
    J0 = cost_function(Q,R,xd,xsol0,u0,tgrid,dt);
    
    phi = compute_phi(A, Q, xsol0, xd, tgrid, Mass);
    gradJ = B.'*phi + R*u0;
    G = innerproduct(gradJ,gradJ,dt);
    
    dx = compute_X(A, 0*X0, B, gradJ, tgrid, Mass);
    H = hessian(Q,R,dx,gradJ,dt);
    step = G/H;

% % This commented part is very useful for debugging. 
% % When the gradient is computed correctly, the two lines should be
% % tangent in step = 0. 
% % When the step size is computed correctly, the two lines should overlap
% % completely. 
%         steps = linspace(0,-step,10);
%         for kk = 1:length(steps)
%             u1 = u0 + steps(kk)*gradJ;
%             xsol1 = compute_X(A, X0, B, u1, tgrid, Mass);
%             J(kk) = cost_function(Q,R,xd,xsol1,u1,tgrid,dt);
%             J2(kk) = cost_function(Q,R,xd,xsol0+steps(kk)*dx,u0+steps(kk)*gradJ,tgrid,dt);
%         end
%         plot(steps, J, steps, J0+G*steps+H/2*steps.^2, steps, J2)
    
    u1 = u0 - step*gradJ;
    xsol1 = xsol0 - step*dx;
    J1 = cost_function(Q,R,xd,xsol1,u1,tgrid,dt);
    
%     if abs(J1-J0) < tol*abs(J0)
    if norm(u1-u0) < tol*norm(u0) %&& norm(xsol1-xsol0) < tol*norm(xsol0)
    norm(u1-u0)/norm(u0)
    ii
    duration = toc;
        return;
    end
    
    u0 = u1;
    %     xsol0 = xsol1;
end

disp('not converged')
duration = toc;
% u0 = [];

end

function J = cost_function(Q,R,xd,x,u,tgrid,dt)
J = 0;
ndt = length(dt);
dx = x - xd(tgrid.');
for ii = 1:ndt
    J = J + dt(ii)*(dx(:,ii+1).'*Q*dx(:,ii+1)/4 + dx(:,ii).'*Q*dx(:,ii)/4 ...
        + u(:,ii).'*R*u(:,ii)/2);
end
end

function out = innerproduct(u,v,dt)
out = 0;
for ii = 1:length(dt)
    out = out + dt(ii)*(u(:,ii).'*v(:,ii));
end
end

function H = hessian(Q,R,dx,du,dt)
H = 0;
ndt = length(dt);
for ii = 1:ndt
    H = H + dt(ii)*(dx(:,ii+1).'*Q*dx(:,ii+1)/2 + dx(:,ii).'*Q*dx(:,ii)/2 ...
        + du(:,ii).'*R*du(:,ii));
end
end