clc;clear all;

p=50;n=50; 
lambda=sort(linspace(10,0.01,p),'descend');
%% Sample covariance
D=diag(sqrt(lambda));  % square root of true spectrum
%generate a rotation matrix to get dependent normals. 
O=orth(randn(p,p));
X=normrnd(0,1,n,p); 
Y=X*D*O;  
sample_spectrum=abs(sort(eig(Y'*Y/n), 'descend'));

% the Quest estimator
addpath('Quest');
Quest_spec= zeros(p,1);
[~,~,tauhat,~,~,~,~, ~,~,~,~]=QuESTimate(Y,0);
Quest_spec=sort(tauhat,'descend');

% ratio
addpath('Opt_Ratio');
ratio_spec = sample_spectrum;
K=10;
ratio_new_all= zeros(p,K);
for i=1:K
    for k=1:K
        [ratio_spec]=RatioL2(ratio_spec,sample_spectrum,n,p,5,20);
        [ratio_spec,~]=Eigen_correction(sample_spectrum,ratio_spec,n,p);
    end
    [ratio_spec]=RatioL2(ratio_spec,sample_spectrum,n,p,1,30);
    ratio_new_all(:,i)= sort(ratio_spec,'descend');
end
ratio_spec = mean(ratio_new_all,2);

t=1;
plot(lambda(t:p), 'b-');hold on;
plot(sample_spectrum(t:p), 'go-');hold on;
% plot(Eigen_spec(t:p), 'r.-.');hold on;
plot(ratio_spec(t:p), 'k.-');hold on;
plot(Quest_spec(t:p), 'r.-.');hold on;

legend('true','sample','Concent','Quest')
%title('phone n=p=50'); % title('linear^{-0.5}');
%saveas(figure(1),'phone_spectrum_50','epsc')
title('linear n=p=50'); % title('linear^{-0.5}');
saveas(figure(1),'Quest_weak','epsc')
