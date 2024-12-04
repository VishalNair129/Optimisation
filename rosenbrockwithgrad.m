function [f,g] = rosenbrockwithgrad(x)
%MATLAB example
%The following code creates the rosenbrockwithgrad function, which includes the gradient as the second output.
% Calculate objective f
f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;

if nargout > 1 % gradient required
    g = [-400*(x(2)-x(1)^2)*x(1) - 2*(1-x(1));
        200*(x(2)-x(1)^2)];
end
end