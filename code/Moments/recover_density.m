function [pdf,t] = recover_density(H,x,f_k)
%   minimize weighted L_1 distance

    k = size(H,2);
    ns = size(x,2);
    V = ones(k,ns);
    V(1,:) = x;
    for i=2:k
        V(i,:) = V(i-1,:).*x;
    end
    A = [V,-eye(k);-V,-eye(k)];
    b = [H';-H'];
    
    Aeq = [ones(1,ns),zeros(1,k)];
    beq= 1;
    f = [zeros(1,ns),1./ (H.*f_k)];
    lb = zeros(ns+k,1);
    %options = optimoptions('linprog','Algorithm','interior-point','Display','iter');
    %x = abs(randn(1,ns));
    %x = x/sum(x);
    [res,fval] = linprog(f,A,b,Aeq,beq,lb,[],[]);
    %fval
    pdf = res(1:ns)';
    t = f*res;
    H'-V*pdf'
end