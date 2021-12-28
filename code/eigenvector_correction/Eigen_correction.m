function [lambda_new, D]=Eigen_correction(lambdasample,Lambdaest_1,n,p)
%% simulated eigenvector 

% 
K=1;
lambda_new_all = zeros(p,K);

for i=1:K
    N=normrnd(0,1,n,p);  
    Y_sim=N*diag(sqrt(Lambdaest_1));  
    w_sim=Y_sim'*Y_sim/n;
    [V0,D0]=eig(w_sim);
    [D,I]=sort(diag(D0),'descend'); 
    V=V0(:,I); 
    lambda_new=diag(V*diag(lambdasample)*V');
    [lambda_new,~]= sort(lambda_new,'descend');
    lambda_new_all(:,i) = abs(lambda_new);
end
lambda_new = mean(lambda_new_all,2);

end