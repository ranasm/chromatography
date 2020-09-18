function [mval] = hill_fun( pp, tt )

t=tt;
mval=ones(size(t));

e=0;  %e>=0
d=1;

a=pp(1);
b=pp(2);
c=pp(3);

idx = find(t > e);
mval(idx)=d-((d-a).*((t(idx)-e).^b))./(c+(t(idx)-e).^b);

end

%Biexponential function (4 parameters)









