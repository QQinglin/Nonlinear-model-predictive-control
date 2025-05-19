% Derivation of Dynamics
clear
syms x1 x2 x3 u1 u2 lambda1 lambda2 lambda3 q1 q2 q3 r1 r2

% inpput vector fields
g1 = [cos(x3);sin(x3);0];
g2 = [0;0;1];

% Dynamics and Costs
f = g1*u1 + g2*u2;
l = q1*x1^4 + q2*x2^2 + q3*x3^4 + r1*u1^4 + r2*u2^4;
Hamiltonian = l + [lambda1,lambda2,lambda3]*f;
dHdx = gradient(Hamiltonian,[x1;x2;x3]);
dHdu = gradient(Hamiltonian,[u1;u2]);


% Parameters for simulation 
param.T = 0.1;                          % can be chosen larger 
param.delta_T = 0.033;
param.N = 21;% ceil(param.T/param.delta_T);  % can be chosen smaller
param.tau = linspace(0,param.T,param.N);
param.t = 30;
param.MPCsteps = ceil(param.t/param.delta_T);
param.fixPointIter = 10;
param.damping = 0.9;

% Parameters for costs (no terminal cost)
param.q1 = 1;
param.q2 = 5;
param.q3 = 0.1;
param.r1 = 0.125*10; % lower values can cause oscillations
param.r2 = 0.0125*10; % lower values can cause oscillations
param.umax = [0.2;1];
param.umin = [-0.2;-1];

% parallel parking 
x0 = [0.75;0;pi];


% Initial guess for state and adjoint trajectory
lambda = ones(param.N,3);  % not zero for parallel parking 
x = zeros(param.N,3);

% MPC-loop
xsim = zeros(param.MPCsteps,3);
usim = zeros(param.MPCsteps,2);
xk = x0;
t = 0;

for k = 1:param.MPCsteps
   
      for j=1:param.fixPointIter
      % solve forward-backward integration 
      xprev = x;
      x = ode2(@(t,x)F(t,x,lambda,param),param.tau,xk);
      % Damping
      x = (1-param.damping)*x + param.damping*xprev;
    
      % solve backward time integration 
      lambdaprev = lambda;
      rev_t = flip(param.tau);
      lambda_r = ode2(@(t,lambda)H(t,x,lambda,param),rev_t,[0;0;0]);
      lambda = flip(lambda_r);
      %Damping
      lambda = (1-param.damping)*lambda + param.damping*lambdaprev;
      end 

      % Apply first input for one sampling step
      u = [h1(x(2,:)',lambda(2,:)',param);h2(x(2,:)',lambda(2,:)',param)];
      x_next = ode2(@(t,x_sim)ffct(t,x_sim,u),[t , t + param.delta_T],xk);
      
      % increment and save simulation state 
      xk = x_next(2,:)';
      xsim(k,:) = x_next(2,:);
      usim(k,:) = u(:,1)';
      t = t + param.delta_T;
end 


% dynamic function 
function xdot = ffct(t,x,u)

xdot = [cos(x(3));
        sin(x(3));
        0] * u(1) + [0;
               0;
               1] * u(2);
end 


% modified Dynamic functions
function xdot = F(t,x,lambda,param)
x3 = x(3);
lambda_d = interp1(param.tau,lambda,t)';
% modified dynamics
xdot = [cos(x3);
        sin(x3);
        0] * h1(x,lambda_d,param) + [0;
               0;
               1] * h2(x,lambda_d,param);
end

function lambdadot = H(t,x,lambda,param)
    q1 = param.q1;
    q2 = param.q2;
    q3 = param.q3;
    
    x_d = interp1(param.tau,x,t)';
    % Modified adjoint dynamics
    lambdadot = - [4 * q1 * x_d(1)^3;
                   2 * q2 * x_d(2);
                   4 * q3 * x_d(3)^3 - lambda(1) * h1(x_d,lambda,param) * sin(x_d(3)) + lambda(2) * h1(x_d,lambda,param) * cos(x_d(3))];    
end 


function u1 = h1(x,lambda,params)
    r1 = params.r1;
    u1max = params.umax(1);
    u1min = params.umin(1);
    u0 = nthroot(- 1/(4*r1) * (lambda(1) * cos(x(3)) + lambda(2) * sin(x(3))),3);

    % projection on input constraints
    if u0 > u1max
        u1 = u1max;
    elseif u0 < u1min
        u1 = u1min;
    else 
        u1 = u0;
    end 
end 

function u2 = h2(~,lambda,params)
    r2 = params.r2;
    u2max = params.umax(2);
    u2min = params.umin(2);
    u0 =nthroot(-1/(4*r2)*(lambda(3)),3);

    % projection on input constraints
    if u0 > u2max
        u2 = u2max;
    elseif u0 < u2min
        u2 = u2min;
    else 
        u2 = u0;
    end 
end 