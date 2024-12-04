A = [1,2,1,0,0; 3,2,0,1,0;1,4,0,0,1];
c = [-3,-4,0,0,0]'; 
x0 = [0,0,12,24,24]'; 
J = [1,2];
I = [3,4,5]; 


xmin = simplex(A,c,x0,I,J)
F_min=c'*xmin

function x = simplex(A,c,x,I,J)
while true
    [m,n] = size(A);

    % (A1) 
    y = A(:,I)'\c(I);

    % (A2) 
    u_J = c(J) - A(:,J)' * y;

    if all(u_J >= 0)
        return
    end

   
    [~,index] = min(u_J);
    r = J(index);

    % (A4) 
    d_I = A(:,I)\A(:,r);

    if all(d_I <= 0)
        return
    end

    % (A5) 
    d_new = d_I;
        for i=1:length(d_new)
            if (d_new(i) < 0)
                d_new(i) = 0;
            end
        end
    [t,s_index] = min(x(I)./d_new);
    s = I(s_index);

    d = zeros(n,1);
    for i = 1:m
        d(I(i)) = d_I(i);
    end

    % (A6) 
    for i=1:n
        if any(I(I~=s) == i) 
            x(i) = x(i) - t*d(i);
        elseif i == r
            x(i) = t;
        else
            x(i) = 0;
        end
    end
    x
    I(s_index) = r;
    J = setdiff(1:n,I);
end
end
