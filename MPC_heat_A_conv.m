clear all
close all
clc

%% Convergence analysis of MPC for the heat equation 
% in which the controller has access to an imperfect plant model with a perturbed A-matrix.  
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

tildeA = A + 0.3*I;

%% define the time grid

That = 12;
nT = ceil(That*100)+1;
tgrid = linspace(0,That,nT);
tgrid2 = tgrid(1:end-1) + diff(tgrid)/2;

U = zeros(1,nT-1);

%% solve the infinite horizon optimal control problem
% tstart = tic;
Pinf = care(full(A),full(B),full(Q));
[Xopt, duration] = compute_X(tildeA-B*(R\(B.'*Pinf)), X0, B, 0*U, tgrid, Mass);
Uopt = -(R\(B.'*Pinf))*Xopt;
Uopt = (Uopt(1:end-1)+Uopt(2:end))/2;
% durationUinf = toc(tstart);

%% Convergence analysis

% set up figures for plotting
figure(5);
set(gcf,'Position',[100 100 800 300])
plot(tgrid2, Uopt, 'Color', 0.3*[1 1 1], 'linewidth', 2)
hold on

figure(6)
hold on

% define the values of tau that should be considered
tau_list = linspace(0.1, 1, 10);
for tauii = 1:length(tau_list)
    
    % set T and tau
    tau = tau_list(tauii);
    T = tau + 4;
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
            [Xoptkk,  duration] = compute_X(tildeA, X0kk, B, Uoptkk(1:(ind3-ind1)), tgrid(ind1:ind3), Mass);
            
            UMPC = [UMPC, Uoptkk(1:(ind3-ind1))];
            ind1 = ind3;
            X0kk = Xoptkk(:,end);
            
            if kk*tau >= 8
                break;
            end
        end
        figure(5)
        plot(tgrid2(1:length(UMPC)), UMPC)
        
        figure(6)
        plot(tgrid2(1:length(UMPC)), UMPC - Uopt(1:length(UMPC)))
        
        % compute the errors
        durationUMPC(tauii) = toc(tstart);
        errors(tauii) = max(abs(UMPC - Uopt(1:length(UMPC))));
    else
        durationUMPC(tauii) = NaN;
        errors(tauii) = NaN;
    end
end

%% display results
lambdainf = eig(A - B*(R\(B'*Pinf)));
muinf = max(real(lambdainf));

figure(7)
ind = 1:length(UMPC);
plot(tgrid2(ind), UMPC(ind), tgrid2(ind), Uopt(ind))

figure(1)
tau_list1 = linspace(0,1,2);
plot(tau_list1, 0.23e-3 + 0.24e-3*tau_list1, 'Color', 0.6*[1 1 1], 'linewidth', 2)
xlabel '\tau'
ylabel '|u_{MPC} - v_\infty|_{L^\infty(0,8)}'
hold on
plot(tau_list, errors, 'k--x', 'linewidth', 2, 'MarkerSize', 10)
legend('0.1 + 0.3 \tau', 'numerical results', 'Location', 'northwest')