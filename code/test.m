clc;clear all;

p=50;n=50; 
flag='linear'; % linear  step  stock  quadratic

if strcmp(flag,'linear') 
    %Use a linear spectrum
    lambda=sort(linspace(10,0.01,p),'descend');
elseif strcmp(flag, 'concave')
    lambda=sort(linspace(10,0.01,p),'descend').^0.3;
elseif strcmp(flag, 'convex')
    lambda=sort(linspace(5,0.01,p),'descend').^2;
elseif strcmp(flag,'step') 
    % use a step spectrum
    lambda=[1+ones(p/2,1);ones(p/2,1)];
    %lambda=[3+ones(p/4,1); 2+ones(p/4,1); 1+ones(p/4,1);ones(p/4,1)]/2;
elseif strcmp(flag, 'sparse')
    lambda=sort(linspace(10,0.01,p),'descend');
    lambda((p/2):p)=0.0001;
elseif strcmp(flag,'stock') 
%      load('DailyReturn800.mat');
%      X0=DailyReturn800(1:n,1:p); %randsample(800,p)
%      w=X0'*X0/n; 
     load('covariance_stock.mat')
     %lambda=sort(repmat(eig(w),1,1),'descend'); 
     lambda=sort(eig(w(1:p,1:p)),'descend');  
elseif strcmp(flag,'amazon') 
     load('Amazon.mat')
     lambda=abs(sort(eig(w(1:p,1:p)),'descend')); 
elseif strcmp(flag,'phone') 
     load('human_phone.mat')
     lambda=abs(sort(eig(w(1:p,1:p)),'descend')); 
elseif strcmp(flag, 'random')
    k=500; breaks=linspace(0,2,k+1); position=linspace(0,1,p);
    coefs = 0.1*randn(k,1); pp = mkpp(breaks,coefs);
     lambda= sort(abs(ppval(pp,position)),'descend'); 
end



%% Sample covariance
D=diag(sqrt(lambda));  % square root of true spectrum
%generate a rotation matrix to get dependent normals. 
O=orth(randn(p,p));
X=normrnd(0,1,n,p); 
Y=X*D*O;  
sample_spectrum=abs(sort(eig(Y'*Y/n), 'descend'));


%% compare

%% the Quest estimator
addpath('Quest');
Quest_spec= zeros(p,1);
[~,~,tauhat,~,~,~,~, ~,~,~,~]=QuESTimate(Y,0);
Quest_spec=sort(tauhat,'descend');

%% moment
% https://weihaokong.github.io/
addpath('Moments');
k = 10;                                        %moments order
for i=1:k                                      %set scaling factor for each moment
    f_k(i) = max(p^(i/2-1),1)/n^(i/2)*(2*i)^(2*i); 
end
Moment_spec=zeros(p,1);
H = compute_moment_by_cycle(Y'*Y,k,p);    %estimate spectral moments
x = 0:(1/p):max(sample_spectrum);                                  %support of recovered measure
[rec_pdf,t] = recover_density(H,x,f_k);        %recover spectral distribution
dis_rec_pdf = pdf2vec(rec_pdf,p);      
% turn pdf into eigenvalues
repeat_time = int8(round(dis_rec_pdf*p));
ind=0;
for i=1:p
    if repeat_time(i)==0
    else
        freq = repeat_time(i); 
        Moment_spec(ind+1:(ind+freq),1)=x(i);
        ind=ind+ freq;
    end
end
Moment_spec=sort(abs(Moment_spec),'descend');


%% ratio
addpath('Opt_Ratio');
if exist('ratio_spec','var') == 0
    ratio_spec = sample_spectrum;
end
ratio_spec = sample_spectrum;
K=10;
ratio_new_all= zeros(p,K);
for i=1:K
    for k=1:K
        [ratio_spec]=RatioL2(ratio_spec,sample_spectrum,n,p,5,20);
        [ratio_spec,~]=Eigen_correction(sample_spectrum,ratio_spec,n,p);
    end
    [ratio_spec]=RatioL2(ratio_spec,sample_spectrum,n,p,1,20);
    ratio_new_all(:,i)= sort(ratio_spec,'descend');
end
ratio_spec = mean(ratio_new_all,2);
if exist('done','var') == 1
    done=done+1
else
    done=1
end


%% plot
t=1;
plot(lambda(t:p), 'b-');hold on;
plot(sample_spectrum(t:p), 'go-');hold on;
% plot(Eigen_spec(t:p), 'r.-.');hold on;
plot(ratio_spec(t:p), 'k.-');hold on;
plot(Quest_spec(t:p), 'g.-.');hold on;
plot(Moment_spec(t:p), 'y.-.');hold on;

legend('true','sample','Concent','Quest','Moment')
%title('phone n=p=50'); % title('linear^{-0.5}');
%saveas(figure(1),'phone_spectrum_50','epsc')
title('phone n=p=50 (remove largest 5)'); % title('linear^{-0.5}');
saveas(figure(1),'phone_spectrum_50_1','epsc')



%% Sample covariance
D=diag(sqrt(ratio_spec));
%generate a rotation matrix to get dependent normals.
for i=1:1
    X=normrnd(0,1,n,p); 
    Y=X*D;  
    sample_spectrum_new=sort(eig(Y'*Y/n), 'descend');
    t=1; 
    plot(lambda(t:p), 'b-');hold on;
    plot(sample_spectrum(t:p), 'go-.');hold on;
    plot(sample_spectrum_new(t:p), 'g.-');hold on;
end



