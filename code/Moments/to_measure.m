function [pos,mas] = to_measure(x)
   x = sort(x); 
   n = size(x,2);
   pos = [];
   mas = [];
   last = -1e10;
   for i=1:n
       if (x(i)==last)
           mas(end)=mas(end)+1/n;
       else
           pos = [pos x(i)];
           mas = [mas 1/n];
           last = x(i);
      end
   end
end