function [mval] = biexp_fun( pp, tt )

t=tt;
mval=ones(size(t));

% e=0;  %e>=0
% d=1;

a=pp(1);
b=pp(2);
c=pp(3);
d=pp(4);
%delay=pp(5);

idx = find(t > d);
% mval(idx)=d-((d-a).*((t(idx)-e).^b))./(c+(t(idx)-e).^b);
% mval(idx)=delay+(a.*exp(1).^((-1).*b.*((-1).*d+t(idx)))+(1+(-1).*a).*exp(1).^(c.*((-1).* ...
%  d+t(idx))));
mval(idx)=(a.*exp(1).^((-1).*b.*((-1).*d+t(idx)))+(1+(-1).*a).*exp(1).^(c.*((-1).* ...
    d+t(idx))));
end