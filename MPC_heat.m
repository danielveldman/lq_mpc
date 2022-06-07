clear all
close all
clc

%% A basic implementation of MPC for the heat equation 
% Daniel Veldman, chair in Dynamics, Control, and Numerics,
% Department of Data Science, FAU Erlangen-Nurnberg. 
% Needs the files:
% 1) compute_control.m
% 2) compute_X.m
% 3) compute_phi.m

% Build the matrix K that is a finite difference discretization of the
% Laplacian. 

L = 1;
N = 11;
x = linspace(0,L,N);
dx = x(2) - x(1);

K = sparse(N,N);
K(1,1) =  1;
K(1,2) = -1;
for ii = 2:N-1
    K(ii,ii-1) = -1;
    K(ii,ii  ) =  2;
    K(ii,ii+1) = -1;
end
K(N,N-1) = -1;
K(N,N  ) =  1;
K = K/dx^2;

% create the matrices in Mass*y'(t) = Ay(t) + Bu(t)
I = speye(N);
A = -K;

B = zeros(N,1);
B(1) = 1;

Mass = I;

% X0 = [erf(2*(x-L/2)).'; zeros(N,1)];
X0 = (x-L/2).'; 

% X0 = [erf(2*(x-L/2)).'; zeros(N,1)];
X0 = [((x-L/2)).'; zeros(N,1)]; % initial condition (position, velocity)

Q = [dx*I*1000, O; O, O]; % weighting matrices and target trajectory
R = 1;
xd =@(t) 0;

%% set the parameters for MPC, the prediction horizon T and control horizon tau

% T   = 4.5;
% tau = 0.5;

tau = 0.5/8;
T   = 4 + tau;

nT = 1000; % number of time points and considered time interval 
That = 4;
tgrid = linspace(0,That,nT); % time grid
tgrid2 = tgrid(1:end-1) + diff(tgrid)/2; % grid of intermediate points (used to plot the controls)

U = zeros(1,nT-1);

%% forward dynamics (for debugging)
% [X, duration] = compute_X(A, X0, B, U, tgrid, Mass);
% figure(1)
% for ii = 1:nT
%     plot(x,X(1:N,ii))
%     title(['t = ', num2str(tgrid(ii))])
%     xlabel 'x'
%     ylabel 'u(t,x)'
%     ylim([-1,1])
%     pause(0.1)
% end

%% adjoint equation (for debugging)
% [phi, duration] = compute_phi(A, Q, X, xd, tgrid, Mass);
% figure(2)
% for ii = 1:nT
%     plot(x,phi(1:N,ii))
%     title(['t = ', num2str(tgrid(ii))])
%     xlabel 'x'
%     ylabel '\phi(t,x)'
%     ylim([-1,1])
%     pause(0.1)
% end

%% classical optimal control (for debugging)
% [Uopt, J0, duration] = compute_control(A, X0, B, U, Q, R, xd, tgrid, Mass);
% [Xopt, duration] = compute_X(A, X0, B, Uopt, tgrid, Mass);
% figure(3)
% plot(tgrid2, Uopt)
% title 'optimal control'
% xlabel 't [s]'
% ylabel 'u_{opt}'

%% compute the infinite horizon optimal control by the solving the ARE
Pinf = care(full(A),full(B),full(Q));
[Xopt, duration] = compute_X(A-B*(R\(B.'*Pinf)), X0, B, 0*U, tgrid, Mass);
Uopt = -(R\(B.'*Pinf))*Xopt;
Uopt = Uopt(1:end-1); %(Uopt(1:end-1)+Uopt(2:end))/2;

% [Uopt, J0, duration] = compute_control(A, X0, B, Uopt, Q, R, xd, tgrid, Mass);
% figure(4)
% for ii = 1:nT
%     plot(x,Xopt(1:N,ii))
%     title(['t = ', num2str(tgrid(ii))])
%     xlabel 'x'
%     ylabel 'u(t,x)'
%     ylim([-1,1])
%     pause(0.1)
% end

%% MPC

%set up the figures in which the control and state will be plotted
figure(5);
plot(tgrid2, Uopt)
hold on
figure(6)
hold on
set(gcf,'Position',[100 100 800 300])
normXopt = zeros(1,size(Xopt,2));
for ii = 1:size(Xopt,2)
    normXopt(ii) = norm(Xopt(:,ii));
end
plot(tgrid, normXopt, 'Color', 0.3*[1 1 1], 'linewidth', 2)

% actual MPC loop. 
UMPC = [];
tstart = 0;
ind1 = 1;
X0kk = X0;
nsteps = ceil(That/tau);
for kk = 1:nsteps
    ind2 = find(tgrid <= (kk-1)*tau + T + eps, 1, 'last'); % find the grid points in the currently considered time interval
    ind3 = find(tgrid <= kk*tau         + eps, 1, 'last');
    tgridkk2 = tgrid(ind1:ind2-1) + diff(tgrid(ind1:ind2));
    
    U0kk = zeros(1,ind2-ind1); 
    % find the optimal control in the current time interval
    [Uoptkk, J0, duration] = compute_control(A, X0kk, B, U0kk, Q, R, xd, tgrid(ind1:ind2), Mass); 
    % compute the state trajectories on [k\tau, (k+1)\tau] and on [k\tau, k\tau + T]
    [Xoptkk,  duration] = compute_X(A, X0kk, B, Uoptkk(1:(ind3-ind1)), tgrid(ind1:ind3), Mass);
    [Xoptkk2, duration] = compute_X(A, X0kk, B, Uoptkk, tgrid(ind1:ind2), Mass);
    
    % plot the obtained controls and state trajectories on the considered
    % subintervals. 
    figure(5)
    plot(tgridkk2, Uoptkk, 'k:')
    if ind3 < ind2
        plot(tgridkk2(1:(ind3-ind1+1)), Uoptkk(1:(ind3-ind1+1)), 'r', 'linewidth', 2)
    else
        plot(tgridkk2(1:(ind3-ind1)), Uoptkk(1:(ind3-ind1)), 'r', 'linewidth', 2)
    end
    
    figure(6)
    normXoptkk2 = zeros(1,size(Xoptkk2,2));
    for ii = 1:length(normXoptkk2)
        normXoptkk2(ii) = norm(Xoptkk2(:,ii));
    end
    normXoptkk  = zeros(1,size(Xoptkk ,2));
    for ii = 1:length(normXoptkk)
        normXoptkk(ii) = norm(Xoptkk(:,ii));
    end
    plot(tgrid(ind1:ind2), normXoptkk2, 'k--')
    plot(tgrid(ind1:ind3), normXoptkk, 'r', 'linewidth', 2)
    
    % update and repeat
    UMPC = [UMPC, Uoptkk(1:(ind3-ind1))];
    ind1 = ind3;
    X0kk = Xoptkk(:,end);
end

%% finalize plotting

fig = figure(7);
plot(tgrid2, UMPC, tgrid2, Uopt)
legend('MPC control', 'Optimal control')
ylabel 'u(t)'
xlabel 't'
saveas(fig, ['MPC_T=', num2str(T*100), '_tau=', num2str(tau*100), '.jpeg']);
% saveas(fig, ['MPC_T=', num2str(T*100), '_tau=', num2str(tau*100), '.fig']);

figure(8)
plot(tgrid2, UMPC - Uopt)
ylabel 'u_{MPC}(t) - u^*_\infty(t)'
xlabel 't'