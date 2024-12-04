MAXit = 1000;
tol1 = 1e-12;
tol2 = 1e-12;
x0 = 1;
%f = @(x) ((x(1)-6)^2 + (x(2)-5)^2) * ((x(1)-1)^2 + (x(2)-1)^2);
%f = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
f=@(x)x^4-3*x^2+2*x;
[sol, h, Fsol,fmin, NTiter,error] = newton(x0, MAXit, tol1, tol2,f);

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


function [sol, h, Fsol,fmin, NTiter,error] = newton(x0, MAXit, tol1, tol2,f)
    d = length(x0);
    v = zeros(d, MAXit);
    v(:, 1) = x0;
    i = 1;

    while i < MAXit
        [F, DF] = CALCfun(f,v(:,i));
        if det(DF)==0
            disp("Hessian not invertible")
            break;
        end
        v(:, i+1) = v(:,i) + -DF\(F);
        h = norm(v(:, i+1) - v(:, i));
        NTiter = i;
        if h < tol1 || norm(F) < tol2
            break;
        end
        i = i + 1;
    end
    
    [Fsol, DF] = CALCfun(f,v(:, i+1));
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



function [F, DF] = CALCfun(f, v)
   num_vars = length(v);
    syms_vars = sym('x', [1 num_vars]); % Create a symbolic vector of variables

    % Create symbolic function from handle
    f_sym = f(syms_vars);

    % Calculate the gradient (F)
    F_sym = gradient(f_sym, syms_vars);
    F_val = subs(F_sym, syms_vars, v'); % Substitute v into the gradient

    % Calculate the Jacobian (DF)
    DF_sym = jacobian(F_sym, syms_vars);
    DF_val = subs(DF_sym, syms_vars, v'); % Substitute v into the Jacobian

    % Convert symbolic expressions to double
    F = double(F_val);
    DF = double(DF_val);

end