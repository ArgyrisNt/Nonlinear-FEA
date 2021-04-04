clear all;
clc;

R = [1;1];
x0 = [4;6];
x = x0;
xall = x;
res = -R;
epsilon = 1e-8;
maxiter = 1e3;
iter = 0;

while any(abs(res)./abs(R) > epsilon)
    J = jacob(x);
    u = - J\res;
    f = Func(x);
    res1 = f-R;
    r0=u'*res;
    r1=u'*res1;
    eta=r0/r1;
    alpha=eta/2;
    x=x+alpha*u;
    f = Func(x);
    res = f-R;
    xall = [xall x];
end
sol = x
plot(xall(1,:))
hold on
plot(xall(2,:))

function F = Func(x)
x1 = x(1);
x2 = x(2);
f1 = x1^2+x2^2-49;
f2 = x1*x2-24;
F = [f1;f2];
end

function J = jacob(x)
J(1,1) = 2*x(1); % df1/dx1
J(1,2) = 2*x(2); % df1/dx2
J(2,1) = x(2); % df2/dx1
J(2,2) = x(1); % df2/dx2
end