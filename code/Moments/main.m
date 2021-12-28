function [dis_rec_pdf,rec_pdf]=main(Y,lambdasample,n,p,k)                                    %main function
     %p dimension
     %n sample size
     %k moments order
     k=32;
    for i=1:k                                      %set scaling factor for each moment
        f_k(i) = max(p^(i/2-1),1)/n^(i/2)*(2*i)^(2*i); 
    end

   H = compute_moment_by_cycle(Y'*Y,k,p);    %estimate spectral moments
    x = linspace(0,max(lambdasample),p);          %support of recovered measure
    [rec_pdf,~] = recover_density(H,x,f_k);        %recover spectral distribution
    dis_rec_pdf = pdf2vec(rec_pdf,p);              %convert to pdf which corresponds to a vector
   

    hold on
    ylim([0 1]);
    stairs([0 x],[0 dis_rec_pdf],'r');
 
    h=legend('Estimated Spectral Distribution','Population Spectral Distribution','Location','southeast');
    hold off
    fprintf('Earthmover distance error: %d\n',mean(dis_rec_err));

end