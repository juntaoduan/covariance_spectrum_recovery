function H = compute_moment_by_cycle(A,k,p)
    n = size(A,1);
    F = triu(A,1);
    F_i = eye(n);
    for i = 1:k
        H(i) = trace(F_i*A)/nchoosek(n,i)/p;
        F_i = F_i*F;
    end    
end