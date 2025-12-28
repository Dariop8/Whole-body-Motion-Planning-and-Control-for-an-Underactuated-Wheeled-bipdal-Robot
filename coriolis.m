function [C, c] = coriolis(M, q, dq)
    n = length(q);
    C = cell(n, 1); 
    for i = 1:n
        Mi = M(:, i); 
        Ci = (1/2) * (jacobian(Mi, q) + jacobian(Mi, q)' - diff(M, q(i)));
        C{i} = Ci;
        
    end
    
    c = sym(zeros(n, 1));
    for i = 1:n
        c(i) = dq' * C{i} * dq;
    end

    c = simplify(c);
end