function [Sigma]=RatioL2(lambdahat,lambdasample,n,p,r,k)
    s=1;Sigma=lambdahat; rate=0.99;
    while s<=r
        D1=diag(sqrt(abs(Sigma)));
        for i=1:k
             Z=normrnd(0,ones(n,p));
             z=D1'*Z'*Z*D1/n;
             N(:,i)=sort(eig(z),'descend')./sort(Sigma,'descend'); %lambdahat; %
        end
        %correction = sum(N,2)./sum(N.^2,2); % sort(  ,'ascend');
        %Sigma0 = lambdasample.*correction;
        Sigma0=sort(lambdasample.*sum(N,2)./sum(N.^2,2),'descend');
        %Sigma0=lambdasample.*sum(N,2)./sum(N.^2,2);
        %Sigma0=sort(lambdasample.*mean(N,2),'descend');
        %Sigma0=lambdasample.*mean(N,2); 
        Sigma=(1-rate)*Sigma + rate*(Sigma0);
    s=s+1;
    end

end