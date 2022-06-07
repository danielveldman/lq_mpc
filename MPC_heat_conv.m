clear all
close all
clc

%% Convergence analysis of MPC for the heat equation 
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

%% define the overall time grid

That = 8;
nT = ceil(That*100)+1;
tgrid = linspace(0,That,nT);
tgrid2 = tgrid(1:end-1) + diff(tgrid)/2;

U = zeros(1,nT-1); % initial guess for the control

%% forward dynamics
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

%% adjoint equation
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

%% classical optimal control
% compute (an very good approximation of) the infinite horizon optimal
% control by solving an optimal control problem on a large time interval
[Uopt, J0, duration] = compute_control(A, X0, B, U, Q, R, xd, tgrid, Mass); 
[Xopt, duration] = compute_X(A, X0, B, Uopt, tgrid, Mass);
% figure(3)
% plot(tgrid2, Uopt2)
% title 'optimal control'
% xlabel 't [s]'
% ylabel 'u_{opt}'

% tstart = tic;
Pinf = care(full(A),full(B),full(Q));
% [Xopt, duration] = compute_X(A-B*(R\(B.'*Pinf)), X0, B, 0*U, tgrid, Mass);
% Uopt = -(R\(B.'*Pinf))*Xopt;
% Uopt = (Uopt(1:end-1)+Uopt(2:end))/2;
% durationUinf = toc(tstart);

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

% choose the sets of values of T and tau that will be considered
T_list   = 0.3 + linspace(0, 10, 11);
tau_list = 0.3; %2*linspace(0.25, 2, 7);

for Tii = 1:length(T_list)
    for tauii = 1:length(tau_list)
        % set T and tau
        T = T_list(Tii);
        tau = tau_list(tauii);
        if tau <= T
            % do MPC
            tstart = tic;
            UMPC = [];
            ind1 = 1;
            X0kk = X0;
            nsteps = ceil(That/tau);
            eps = 1e-10;
            for kk = 1:nsteps
                ind2 = find(tgrid <= (kk-1)*tau + T + eps, 1, 'last');
                ind3 = find(tgrid <= kk*tau         + eps, 1, 'last');
                tgridkk2 = tgrid(ind1:ind2-1) + diff(tgrid(ind1:ind2));
                
                U0kk = zeros(1,ind2-ind1);
                [Uoptkk, J0, duration] = compute_control(A, X0kk, B, U0kk, Q, R, xd, tgrid(ind1:ind2), Mass);
                [Xoptkk,  duration] = compute_X(A, X0kk, B, Uoptkk(1:(ind3-ind1)), tgrid(ind1:ind3), Mass);
                
                UMPC = [UMPC, Uoptkk(1:(ind3-ind1))];
                ind1 = ind3;
                X0kk = Xoptkk(:,end);
            end
            % save errors and computational time
            durationUMPC(Tii,tauii) = toc(tstart);
            errors(Tii,tauii) = max(abs(UMPC - Uopt(1:length(UMPC))));
        else
            durationUMPC(Tii,tauii) = NaN;
            errors(Tii,tauii) = NaN;
        end
    end
end

%% display results
lambdainf = eig(A - B*(R\(B'*Pinf)));
muinf = max(real(lambdainf));

figure(1)
T_list1 = linspace(0,10,100);
plot(T_list1, 0.002*exp(2*muinf*T_list1), 'Color', 0.6*[1 1 1], 'linewidth', 2)
xlabel 'T'
ylabel '|u_{MPC} - u^*_\infty|_{L^\infty(0,8)}'
hold on
plot(T_list-tau_list, errors, 'k--x', 'linewidth', 2, 'MarkerSize', 10)
legend('0.002 exp(-2 \mu_\infty (T-\tau))', 'numerical results')