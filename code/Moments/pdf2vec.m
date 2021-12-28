function dis_pdf = pdf2vec(pdf,p)
    t = size(pdf,2);
    dis_pdf = zeros(1,t);
    m=0;
    for i=1:t
        m = m+pdf(i);
        while(m>=1/p)
            m = m-1/p;
            dis_pdf(i)=dis_pdf(i)+1/p;
        end
    end
end