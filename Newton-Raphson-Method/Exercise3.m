% Exercise 3
%
% Newton-Raphson Method for solving the system of nonlinear equations with
% incremental load
%               x1^2+x2^2-P1(t)=0
%               x1*x2-P2(t)=0
% where 0<=P1(t)<=50 and 0<=P2(t)<=25


clear all;
clc;

x0 = [0.01;0]; % Initial value
P = [0;0]; % Initial load
P1_all = P(1);
sol(:,1) = Func(x0,P); % Evaluate function at initial value x0 and load P
nvar = length(sol);
h = 1e-4 .* ones(nvar, 1); % Step for the Newton-Raphson Method
Pmax = 50; % Set the maximum value of load
% incr = 10;
incr = 100;
% incr = 1000;

for k = 1:incr
    P(1) = k*Pmax/incr; % Update load in x-direction
    P(2) = k*Pmax/(2*incr);% Update load in y-direction
    iter = 0; % Initialize iteration counter
    maxiter = 1e3; % Set maximum number of iterations
    error1 = 1e8; % Set errors
    error = zeros(maxiter, 1);
    x = x0;
    
    % Run until error1 is not longer greater than 1e-12 or maxiter has been
    % reached.
    while (error1 > 1e-12 & iter < maxiter)
        iter = iter+1;
        f = Func(x,P); % Evaluate function at x=(x1,x2)
        J = jacob(x); % Evaluate Jacobian at x=(x1,x2)  
        x = x - J\f;  % Update x=(x1,x2) 
        error1 = sqrt(sum((J\f).^2)); % Update errors
        error(iter) = sqrt(sum(f.^2));
    end
    % Print result and errors
    P1_all = [P1_all P(1)]; % Store all values of P1 component of load P
    sol(:,k+1) = x; % Store the solutions for different amount of load P
    error = error(1:iter);
end
% Draw the x1-P1 curve
plot(sol(1,:),P1_all)

% Print results
A = [P1_all;sol(1,:)];
fprintf('%6s %12s\r\n','Px','x');
fprintf('%6.2f %12.8f\r\n',A);

function F = Func(x,P)
x1 = x(1);
x2 = x(2);
P1 = P(1);
P2 = P(2);
f1 = x1^2+x2^2-P1;
f2 = x1*x2-P2;
F = [f1;f2];
end

function J = jacob(x)
J(1,1) = 2*x(1); % df1/dx1
J(1,2) = 2*x(2); % df1/dx2
J(2,1) = x(2); % df2/dx1
J(2,2) = x(1); % df2/dx2
end