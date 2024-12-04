MAXit = 1000;
tol1 = 1e-12;
tol2 = 1e-12;
x0 = -1;

a = 1;
b = 100;
%f= @(x) (a - x(1))^2 + b*(x(2) - x(1)^2)^2;
%f=@(x) x(1)^2+2*x(1)*x(2)+x(2)^2;
%f = @(x) ((x(1)-6)^2 + (x(2)-5)^2) * ((x(1)-1)^2 + (x(2)-1)^2);
f=@(x)x^4-3*x^2+2*x;

tau=1e-12;

[sol, H, Fsol,fmin, NTiter,error] = newton(x0, MAXit, tol1, tol2,f,tau);

% Display the solution
disp(['Solution: ', num2str(sol')]);
disp(['Function value at solution: ', num2str(Fsol')]);
disp(['Number of iterations: ', num2str(NTiter)]);
disp(['Minima of the function: ', num2str(fmin)]);


figure(1)
plot(1:NTiter,error,'-o')
xlabel("No of Iterations")
ylabel("Error")
title("Error vs Iterations")
grid on

figure(2)
plot(1:NTiter,log(error),'-o')
xlabel("No of Iterations")
ylabel("Error (log Scale)")
title("Error vs Iterations")
grid on


function [sol, H, Fsol,fmin, NTiter,error] = newton(x0, MAXit, tol1, tol2,f,tau)
    d = length(x0);
    v = zeros(d, MAXit);
    v(:, 1) = x0;
    i = 1;
    while i < MAXit
        [F, DF] = CALCfun(f,v(:,i),tau);
        if det(DF)==0
            disp("Hessian not invertible")
            break;
        end
        v(:, i+1) = v(:,i) + -DF\(F);
        H = norm(v(:, i+1) - v(:, i));
        NTiter = i;
        if H < tol1 || norm(F) < tol2
            break;
        end
        i = i + 1;
    end
    [Fsol, DF] = CALCfun(f,v(:,i+1),tau);
    deT=det(DF);
    if deT > 0 && DF(1,1) > 0
        disp('Minima attained');
    else
        if deT > 0 && DF(1,1) < 0
            disp('We got a Local maximum, Change Initial guess');
        else
            if deT < 0
                disp('We got a Saddle point, Change Initial guess');
            end
        end
    end
    v = v(:, 1:i+1);
    sol = v(:, end); 
    fmin = f(sol); % Calculate the function value at the solution
    error=zeros(1,length(v)-1);
    for i=1:length(v)-1
    error(i)=norm(v(:, i+1) - v(:, i));
    end
    

end


function [F, DF] = CALCfun(f,v,tau)
   num_vars = length(v);
   syms_vars = sym('x', [1 num_vars]); % Create a symbolic vector of variables

    % Create symbolic function from handle
    f_sym = f(syms_vars);
    % Calculate the gradient (F)
    F_sym = gradient(f_sym, syms_vars);
    F_val = double(subs(F_sym, syms_vars, v')); % Substitute v into the gradient and convert to numeric
    F = F_val;
    
    DF = zeros(num_vars);  % Initialize the Jacobian matrix

    for k = 1:num_vars
        for l = 1:num_vars
            % Perturb the l-th variable
            z_temp = v;
            z_temp(l) = z_temp(l) + tau;
            F_temp = double(subs(F_sym, syms_vars, z_temp')); % Gradient at z_temp
            % Secant method for Jacobian: difference of gradient components
            DF(k, l) = (F_temp(k) - F(k)) / tau; % Secant method approximation
        end
    end 
end