beta=rand(); 
epsilon=10^(-12); 
Jmax=1000; 
X0=[0.000001,0.00001];
fmin=0;
f=@(x,y) (x-y)^2;
grad_f=@(x,y) [2*(x-y);2*(y-x)];
[X_min,f_min]=Proj_gradient_descent(f,fmin,grad_f,beta,epsilon,Jmax,X0);
X_min
f_min
function [X_min,f_min]=Proj_gradient_descent(f,fmin,grad_f,beta,epsilon,Jmax,X0)
k=0;
norm_f=abs(f(X0(1),X0(2))-fmin);
while norm_f>epsilon && k <Jmax
s=0;
grad=grad_f(X0(1),X0(2));
p=Proj_circ(X0(1)-(beta^s)*grad(1),X0(2)-(beta^s)*grad(2));
while f(X0(1),X0(2))<f(p(1),p(2))
s=s+1;
p=Proj_circ(X0(1)-(beta^s)*grad(1),X0(2)-(beta^s)*grad(2));
end
X0=X0-beta^s*grad_f(X0(1),X0(2));
X0=Proj_circ(X0(1),X0(2));
norm_f=norm(f(X0(1),X0(2))-fmin);
k=k+1;
end
X_min=X0;
f_min=f(X0(1), X0(2));
end

function p=Proj_circ(x,y)
X=[x;y];
zero=[0;0];
if X==zero
p=[1;0];
else
p=[x/sqrt(x^2+y^2);y/sqrt(x^2+y^2)];
end
end

