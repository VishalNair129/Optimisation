
broyden_method_with_plots();



function broyden_method_with_plots()
    
    F = @(x) [ (x(1) + 3)*(x(2)^2 - 7) + 18; sin(x(2) * exp(x(1)) - 1) ];

    % Initial guess
    x0 = [1; 1];  % Starting point

    % Tolerance and maximum number of iterations
    tol = 1e-6;
    max_iter = 500;
    beta = 0.5;  % Parameter for line search
    
    % Run the Broyden method
    [x_sol, num_iterations, x_vals, F_vals] = broyden_method(F, x0, tol, max_iter, beta);
    [x_sol2, num_iterations2, x_vals2, F_vals2] = broyden_method(F, x0, tol, max_iter, beta);
    % Plot the solution (xk, yk) at each iteration
    figure;
    hold on
    plot(x_vals(1,:), x_vals(2,:), 'o-');
    plot(x_vals2(1,:), x_vals2(2,:), 'g-');
    xlabel('x');
    ylabel('y');
    legend('Broyden method','Broyden method with SMW')
    title('Solution (x_k, y_k) at each iteration');

    % Plot the norm F(xk, yk) at each iteration
    figure;
    semilogy(1:num_iterations, F_vals, 'o-');
    hold on
    semilogy(1:num_iterations2,F_vals2,'g-')
    xlabel('Iteration');
    ylabel('Norm of F(x_k, y_k)');
    title('Norm of F(x_k, y_k) at each iteration');

    % Runtime comparison
    % Direct matrix inversion
    tic;
    broyden_method(F, x0, tol, max_iter, beta);
    time_direct = toc;

    % Sherman-Morrison-Woodbury
    tic;
    broyden_method_smw(F, x0, tol, max_iter, beta);
    time_smw = toc;

    fprintf('Runtime using direct matrix inversion: %.4f seconds\n', time_direct);
    fprintf('Runtime using Sherman-Morrison-Woodbury: %.4f seconds\n', time_smw);
    fprintf('Solution given by Broyden (x_k, y_k) = (%.6f, %.6f)\n', x_sol(1,end), x_sol(2,end));
    fprintf('Solution given by Broyden SMW (x_k, y_k) = (%.6f, %.6f)\n', x_sol2(1,end), x_sol2(2,end));
end

function [x, k, x_vals, F_vals] = broyden_method(F, x0, tol, max_iter, beta)
    % Initialize
    x = x0;
    n = length(x);
    H = eye(n);  % Initial inverse Jacobian approximation
    x_vals = x;           % Store x values for plotting
    F_vals = norm(F(x));  % Store F values for plotting

    for k = 1:max_iter
        Fx = F(x);
        if norm(Fx) < tol
            return;
        end
        
        delta_x = -H \ Fx;  % Solve for the step

        % Line search to choose the step size
        t = 1;
        while norm(F(x + t * delta_x)) >= norm(Fx) && t > tol
            t = beta * t;
        end

        % Update the estimate
        x_new = x + t * delta_x;
        
        s = x_new - x;
        if norm(s) < tol  % Stop condition to avoid division by zero
            return;
        end

        y = F(x_new) - Fx;
        
        % Update the inverse Jacobian approximation
        H = H + ((y - H * s) * s') / (s' * s);
              
        x = x_new;

        % Store values for plotting
        x_vals(:, end+1) = x;
        F_vals(end+1) = norm(F(x));
    end

    warning('Maximum number of iterations reached');
end

function [x, k] = broyden_method_smw(F, x0, tol, max_iter, beta)
    % Initialize
    x = x0;
    n = length(x);
    H = eye(n);  % Initial inverse Jacobian approximation

    for k = 1:max_iter
        Fx = F(x);
        if norm(Fx) < tol
            return;
        end
        
        delta_x = -H * Fx;  % Solve for the step

        % Line search to choose the step size
        t = 1;
        while norm(F(x + t * delta_x)) >= norm(Fx) && t > tol
            t = beta * t;
        end

        % Update the estimate
        x_new = x + t * delta_x;
        
        s = x_new - x;
        if norm(s) < tol  % Stop condition to avoid division by zero
            return;
        end

        y = F(x_new) - Fx;
        
       % Update the inverse Jacobian approximation using Sherman-Morrison-Woodbury
 
        Hy = H \ y;
        Hs = H \ s;
        H = H\eye(n) + (s - Hy) * Hs' / (s' * Hs + 1);
        
        % Update the current estimate
        x = x_new;
    end

    warning('Maximum number of iterations reached');
end

