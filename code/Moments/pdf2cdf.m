function cdf = pdf2cdf(pdf)
    last = 0;
    for i = 1:size(pdf,2)
       last = last + pdf(i);
       cdf(i) = last;
    end
    
end