gradient_descent_comparison();

function gradient_descent_comparison()
    x0 = 0.0;
    y0 = 0.0;
    alpha = 0.01;
    tol = 1e-6;
    max_iter = 10000;

    [xy_const, iter_const] = gradient_descent_constant_step([x0; y0], alpha, tol, max_iter);
    fprintf('Constant Step Size: x = %.6f, y = %.6f, iterations = %d\n', xy_const(1), xy_const(2), iter_const);

    [xy_wolfe, iter_wolfe] = gradient_descent_wolfe([x0; y0], tol, max_iter);
    fprintf('Wolfe Step Size Rule: x = %.6f, y = %.6f, iterations = %d\n', xy_wolfe(1), xy_wolfe(2), iter_wolfe);
end

function f = cost_function(xy)
    x = xy(1);
    y = xy(2);
    f = (3*x - exp(x)*y - 2)^2 + (x + y - 1)^2;
end

function grad = gradient(xy)
    x = xy(1);
    y = xy(2);
    df_dx = 2 * (3*x - exp(x)*y - 2) * (3 - y*exp(x)) + 2 * (x + y - 1);
    df_dy = 2 * (3*x - exp(x)*y - 2) * (-exp(x)) + 2 * (x + y - 1);
    grad = [df_dx; df_dy];
end

function [xy, iter] = gradient_descent_constant_step(xy0, alpha, tol, max_iter)
    xy = xy0;
    for iter = 1:max_iter
        grad = gradient(xy);
        if norm(grad) < tol
            break;
        end
        xy = xy - alpha * grad;
    end
end

function [xy, iter] = gradient_descent_wolfe(xy0, tol, max_iter)
    theta = 0.9;
    xy = xy0;
    for iter = 1:max_iter
        grad = gradient(xy);
        if norm(grad) < tol
            break;
        end
        pk = -grad;
        alpha = 1.0; % initial step size
        while ~wolfe_condition(@gradient, xy, pk, alpha, grad, theta)
            alpha = alpha * 0.5; % backtracking
        end
        xy = xy + alpha * pk;
    end
end


function result = wolfe_condition(grad_f, xy, pk, alpha, grad, theta)
    grad_new = grad_f(xy + alpha * pk);
    result = abs(pk' * grad_new) <= theta * abs(pk' * grad);
end
