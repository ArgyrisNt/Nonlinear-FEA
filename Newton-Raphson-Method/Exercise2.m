% Exercise 2
%
% Newton-Raphson Method for solving the system of nonlinear 
% equations:
%               x1^2+x2^2-50=0
%               x1*x2-25=0
% where x1,x2>0


clear all;
clc;

error1 = 1e8; % Initialize error for the while loop
x0 = [1;0]; % Initial value
x = x0;
iter = 0; % Initialize iteration counter
maxiter = 1e3; % Set maximum number of iterations
error = zeros(maxiter, 1);
sol = Func(x0); % Evaluate function an initial value x0
nvar = length(sol);
h = 1e-4 .* ones(nvar, 1); % Step for the Newton-Raphson Method
xall = x;

while error1 > 1e-12 % Run until error1 is not longer greater than 1e-12
    iter = iter+1;
    f = Func(x); % Evaluate function at x=(x1,x2)
    J = jacob(x); % Evaluate Jacobian at x=(x1,x2)  
    x = x - J\f;  % Update x=(x1,x2) 
    error1 = sqrt(sum((J\f).^2)); % Update errors
    error(iter) = sqrt(sum(f.^2));
    if iter == maxiter % Check if maxiter has been reached
        warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
        break;
    end
    xall = [xall x]; % Store all values of x
end
% Print result and errors
sol = x;
error = error(1:iter);

% Show that solution converges
plot(xall(1,:))
hold on
plot(xall(2,:))

% Print results
A = [xall(1,:);xall(2,:)];
fprintf('%6s %12s\r\n','x','y');
fprintf('%6.2f %12.8f\r\n',A);

function F = Func(x)
x1 = x(1);
x2 = x(2);
f1 = x1^2+x2^2-50;
f2 = x1*x2-25;
F = [f1;f2];
end

function J = jacob(x)
J(1,1) = 2*x(1); % df1/dx1
J(1,2) = 2*x(2); % df1/dx2
J(2,1) = x(2); % df2/dx1
J(2,2) = x(1); % df2/dx2
end