subgradient_algorithm();
function subgradient_algorithm()
    % Parameters
    x0 = [1.0; 1.0];  % Initial point
    epsilon = 1e-6;  % Tolerance
    beta = 0.9999;  % Step size parameter for beta^k
    beta2 = 0.9999; % Step size parameter for beta/(k+1)
    max_iter = 1000000;  % Maximum number of iterations

    % Run the algorithm with beta^k
    fprintf('Subgradient Algorithm with beta^k:\n');
    [optimal_x, optimal_f] = subgradient_algorithm_core(x0, epsilon, beta, max_iter);
    fprintf('Optimal x: (%.6f, %.6f)\n', optimal_x(1), optimal_x(2));
    fprintf('Optimal f(x): %.6f\n\n', optimal_f);
    
    % Run the algorithm with beta/(k+1)
    fprintf('Subgradient Algorithm with beta/(k+1):\n');
    [optimal_x, optimal_f] = subgradient_algorithm_core2(x0, epsilon, beta2, max_iter);
    fprintf('Optimal x: (%.6f, %.6f)\n', optimal_x(1), optimal_x(2));
    fprintf('Optimal f(x): %.6f\n', optimal_f);
end

function [xk, fk] = subgradient_algorithm_core(x0, epsilon, beta, max_iter)
    xk = x0;
    k = 0;
    
    while norm(f(xk)) > epsilon && k < max_iter
        sk = subgradient(xk);
        dk = sk / norm(sk);
        tk = beta^k;
        xk = xk - tk * dk;
        k = k + 1;
       
    end
    
    fk = f(xk);
    k
end

function [xk, fk] = subgradient_algorithm_core2(x0, epsilon, beta, max_iter)
    xk = x0;
    k = 0;
    
    while norm(f(xk)) > epsilon && k < max_iter
        sk = subgradient(xk);
        dk = sk / norm(sk);
        tk = beta / (k + 1);
        xk = xk - tk * dk;
        k = k + 1;
    end
    
    fk = f(xk);
    k
end

function val = f(x)
    x1 = x(1);
    x2 = x(2);
    f1 = x1^2 + (x2 - 1)^2 + x2 - 1;
    f2 = -x1^2 - (x2 - 1)^2 + x2 + 1;
    val = max(f1, f2);
end

function sg = subgradient(x)
    x1 = x(1);
    x2 = x(2);
    f1 = x1^2 + (x2 - 1)^2 + x2 - 1;
    f2 = -x1^2 - (x2 - 1)^2 + x2 + 1;
    if f1 > f2
        sg = [2*x1; 2*(x2 - 1) + 1];
    else
        sg = [-2*x1; -2*(x2 - 1) + 1];
    end
end
