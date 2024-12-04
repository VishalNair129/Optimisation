
f = @(x) x(1)^2 + 3*x(2)^2;
h = @(x) 2*x(1)^2 + x(2)^2 - 8;
grad_f = @(x) [2*x(1); 6*x(2)];
grad_h = @(x) [4*x(1); 2*x(2)];

F = @(x, mu) [grad_f(x) + mu * grad_h(x); h(x)];

jac_F = @(x, mu) [
    2,    0,  4*x(1);
    0,    6,  2*x(2);
    4*x(1), 2*x(2), 0
];


x0 = [0; 1]; % Starting point for x
mu0 = 1;     % Starting point for mu
tol = 1e-6;  % Tolerance 

[x_min, mu_min, num_iterations] = newton_method_multi(x0, mu0, tol,F,jac_F);

% Display results
disp('Found minimum for x:');
disp(x_min);
disp('Found minimum for mu:');
disp(mu_min);
disp('Number of iterations:');
disp(num_iterations);

function [x, mu, k] = newton_method_multi(x0, mu0, tol,F,jac_F)

    x = x0;
    mu = mu0;
    k = 0;

    while norm(F(x, mu)) >= tol
        if all(F(x, mu) == 0)
            error('F(x, mu) is zero. Aborting.');
        end

        % Solve for the Newton direction
        delta = -jac_F(x, mu) \ F(x, mu);
        delta_x = delta(1:2);
        delta_mu = delta(3);

        % Update step
        x = x + delta_x;
        mu = mu + delta_mu;

        % Increment iteration counter
        k = k + 1;
    end
end