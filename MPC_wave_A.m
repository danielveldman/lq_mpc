clear all
close all
clc

%% A basic implementation of MPC for the wave equation 
% in which the controller has access to an imperfect plant model with a perturbed A-matrix.  
% this is a minor modification of MPC_wave.m
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

% Build the matrix A that represents x''(t) = Kx(t) 
% as a first order system y'(t) = A y(t), 
% where y(t) = [x(t), x'(t)]. 
O = 0*K;
I = speye(size(K));
A = [O, I; -K, O];  % A-matrix used by the MPC controller

B = zeros(2*N,1); % create the matrices in Mass*y'(t) = Ay(t) + Bu(t)
B(N+1) = 1;       % (the control is applied at the last node)

Mass = speye(2*N); 

% X0 = [erf(2*(x-L/2)).'; zeros(N,1)];
X0 = [((x-L/2)).'; zeros(N,1)]; % initial condition (position, velocity)

Q = [dx*I*1000, O; O, O]; % weighting matrices and target trajectory
R = 1;
xd =@(t) 0;

tildeA = A + [O, O; O, 0.3*I];  % A-matrix of the system that is being controlled. 

%% set the parameters for MPC, the prediction horizon T and control horizon tau

% T   = 4.5;
% tau = 0.5;

tau = 0.5/8;
T   = 4 + tau;

That = 12; % length of the time interval considered in the whole sim
nT = ceil(That*100)+1;
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
% [Uopt2, J0, durationUopt] = compute_control(A, X0, B, U, Q, R, xd, tgrid, Mass);
% figure(3)
% plot(tgrid2, Uopt2)
% title 'optimal control'
% xlabel 't [s]'
% ylabel 'u_{opt}'

%% compute the infinite horizon optimal control and the MPC limit by solving the AREs

tstart = tic;
Pinf = care(full(A),full(B),full(Q));
Acl = A-B*(R\(B.'*Pinf));
lambdainf = eig(Acl);
muinf = max(real(lambdainf));

Qinf = care(full(tildeA),full(B),full(Q));
tildeAcl = A-B*(R\(B.'*Qinf));

[Zopt, duration] = compute_X(tildeA-B*(R\(B.'*Pinf)), X0, B, 0*U, tgrid, Mass);
Vopt = -(R\(B.'*Pinf))*Zopt;
Vopt = (Vopt(1:end-1)+Vopt(2:end))/2;

[Xopt, duration] = compute_X(tildeAcl, X0, B, 0*U, tgrid, Mass);
Uopt = -(R\(B.'*Pinf))*Xopt;
Uopt = (Uopt(1:end-1)+Uopt(2:end))/2;

durationUinf = toc(tstart);

% [Uopt3, J0, durationUopt] = compute_control(A, X0, B, Uopt, Q, R, xd, tgrid, Mass);
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

% set up figures for plotting the control and state trajectories
UMPC = [];
tstart = 0;
ind1 = 1;
X0kk = X0;
nsteps = ceil(That/tau);
figure(5);
set(gcf,'Position',[100 100 800 300])
hU = plot(tgrid2, Uopt, 'Color', 0.3*[1 1 1], 'linewidth', 2);
hold on
hV = plot(tgrid2, Vopt, 'Color', 0.6*[1 1 1], 'linewidth', 2);
% plot(tgrid2, Uopt2, 'Color', 0.3*[0 1 1], 'linewidth', 2)
% plot(tgrid2, Uopt3, 'Color', 0.3*[0 0 1], 'linewidth', 2)

figure(6)
hold on
set(gcf,'Position',[100 100 800 300])
normZopt = zeros(1,size(Zopt,2));
normXopt = zeros(1,size(Xopt,2));
for ii = 1:size(Xopt,2)
    normZopt(ii) = norm(Zopt(:,ii));
    normXopt(ii) = norm(Xopt(:,ii));
end
hX = plot(tgrid, normXopt, 'Color', 0.3*[1 1 1], 'linewidth', 2);
hZ = plot(tgrid, normZopt, 'Color', 0.6*[1 1 1], 'linewidth', 2);

% actual MPC loop. 
tstart = tic;
eps = 1e-10;
for kk = 1:nsteps    
    ind2 = find(tgrid <= (kk-1)*tau + T + eps, 1, 'last'); % find the grid points in the currently considered time interval
    ind3 = find(tgrid <= kk*tau         + eps, 1, 'last');
    tgridkk2 = tgrid(ind1:ind2-1) + diff(tgrid(ind1:ind2));
    ind4 = min([find(tgridkk2 <= kk*tau + eps, 1, 'last') + 1, length(tgridkk2)]);
    
    U0kk = zeros(1,ind2-ind1);
    % find the optimal control in the current time interval (based on the matrix A in the system model)
    [Uoptkk, J0, duration] = compute_control(A, X0kk, B, U0kk, Q, R, xd, tgrid(ind1:ind2), Mass);
    % compute the state trajectories on [k\tau, (k+1)\tau] and on [k\tau,
    % k\tau + T] (based on the true system matrix tildeA)
    [Xoptkk,  duration] = compute_X(tildeA, X0kk, B, Uoptkk(1:(ind3-ind1)), tgrid(ind1:ind3), Mass);
    [Xoptkk2, duration] = compute_X(tildeA, X0kk, B, Uoptkk, tgrid(ind1:ind2), Mass);
    
    % plot the obtained controls and state trajectories on the considered
    % subintervals. 
    figure(5)
%     hUp = plot(tgridkk2, Uoptkk, 'k--');
    hUMPC = plot(tgridkk2(1:ind4), Uoptkk(1:ind4), 'r', 'linewidth', 2);
    
    figure(6)
    normXoptkk2 = zeros(1,size(Xoptkk2,2));
    for ii = 1:length(normXoptkk2)
        normXoptkk2(ii) = norm(Xoptkk2(:,ii));
    end
    normXoptkk  = zeros(1,size(Xoptkk ,2));
    for ii = 1:length(normXoptkk)
        normXoptkk(ii) = norm(Xoptkk(:,ii));
    end
%     hXp = plot(tgrid(ind1:ind2), normXoptkk2, 'k--');
    hXMPC = plot(tgrid(ind1:ind3), normXoptkk, 'r', 'linewidth', 2);
    
    % stop if a sufficiently large time interval has been considered
    if kk*tau >= 7
        break
    end
    
    % update and repeat
    ind1 = ind3;
    UMPC = [UMPC, Uoptkk(1:ind4)];
    X0kk = Xoptkk(:,end);
end
durationUMPC = toc(tstart);

%% finalize plotting

fig = figure(5);
% plot(tgrid2, UMPC, tgrid2, Uopt)
set(gca,'FontSize',18) 
xlim([0 7])
ylim([-10 15])
grid on
lgd = legend([hU, hV, hUMPC], {'Infinite-horizon control u_\infty^*', 'MPC limit v_\infty', 'MPC control u_{MPC}'});
lgd.Position = lgd.Position + [0 0.1 0 0];
ylabel 'u(t)'
xlabel 't'
name = ['MPCA_T=', num2str(T*10000), '_tau=', num2str(tau*10000), '.jpeg'];
print('-djpeg','-r500', name)
% saveas(fig, ['MPC_T=', num2str(T*10000), '_tau=', num2str(tau*10000), '.fig']);

fig = figure(6);
% plot(tgrid2, UMPC, tgrid2, Uopt)
xlim([0 7])
ylim([0 5])
grid on
% legend([hX, hZ, hXMPC, hXp], {'Infinite-horizon trajectory y_\infty^*', 'MPC limit z_\infty', 'MPC trajectory y_{MPC}', 'MPC trajectories on [k\tau, k\tau+T]'})
legend([hX, hZ, hXMPC], {'Infinite-horizon trajectory y_\infty^*', 'MPC limit z_\infty', 'MPC trajectory y_{MPC}'})
ylabel '|y(t)|'
xlabel 't'
name = ['MPCAX_T=', num2str(T*10000), '_tau=', num2str(tau*10000), '.jpeg'];
print('-djpeg','-r500', name)
% saveas(fig, ['MPC_T=', num2str(T*10000), '_tau=', num2str(tau*10000), '.fig']);