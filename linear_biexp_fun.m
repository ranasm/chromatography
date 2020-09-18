function [mval] = linear_biexp_fun( pp, tt )

t=tt;
mval=ones(size(t));

A1=pp(1);
A2=1;
A3=pp(2);
lambda1=pp(3);
lambda2=pp(4);
tau=pp(5);
pre= t<tau;
post= t>=tau;

mval(pre)=A1*t(pre)+A2;
mval(post)=(A1*tau+A2+A3)/2*exp(-lambda1*(t(post)-tau))+(A1*tau+A2-A3)/2*exp(-lambda2*(t(post)-tau)) ;

end