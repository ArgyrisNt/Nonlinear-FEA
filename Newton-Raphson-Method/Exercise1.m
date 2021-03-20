% Exercise 1
%
% Newton-Raphson Method for solving a nonlinear equation


clear all;
clc;

x0 = 0.1; % Initial value
P = 0; % Initial load
sol = Func(x0,P); % Evaluate function at initial value and load P
nvar = length(sol);
h = 1e-4 .* ones(nvar, 1); % Step for the Newton-Raphson Method
Pmax = 1000; % Set the maximum value of load
incr = 100;
P_all = P;

for k = 1:incr
    P = k*Pmax/incr; % Update load   
    iter = 0; % Initialize iteration counter
    maxiter = 1e3; % Set maximum number of iterations
    error1 = 1e8; % Set errors
    error = zeros(maxiter, 1);
    x = x0;
    % Run while loop until error1 is not longer greater than 1e-12 or
    % maxiter has been reached.
    while (error1 > 1e-12 & iter < maxiter)
        iter = iter+1;
        f = Func(x,P); % Evaluate function at x
        df = derf(x); % Evaluate derivative of function at x
        x = x - f/df;  % Update x
        error1 = sqrt(sum((f/df).^2)); % Update errors
        error(iter) = sqrt(sum(f.^2));
    end
    % Print result and errors
    sol(k+1) = x; % Store the solutions for different amount of load P
    error = error(1:iter);
    P_all = [P_all P]; % Store all values of load P
end
% Draw the x-P curve
plot(sol,P_all)

% Print results
A = [P_all;sol];
fprintf('%6s %12s\r\n','P','x');
fprintf('%6.2f %12.8f\r\n',A);

function F = Func(x,P)
F = (-1+1/sqrt(1-2*x*sin(15)+x^2))*(sin(15)-x)-P;
end

function df = derf(x)
df = (-1/2*(1-2*x*sin(15)+x^2)^(-3/2)*(2*x-2*sin(15)))*(sin(15)-x)-(-1+1/sqrt(1-2*x*sin(15)+x^2));
end
