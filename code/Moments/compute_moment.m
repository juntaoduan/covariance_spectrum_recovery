function beta = compute_moment(D,k)
    n = size(D,2);
    beta = [];
    for i=1:k
        beta = [beta sum(D.^i)/n];
    end
end