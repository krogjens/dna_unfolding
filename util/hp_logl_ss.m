function logl = hp_logl_ss(data,theta,n)

gamma = theta(1);
xc = data.xc;
L = data.L;
var = 2/(gamma * L)*n;

logl = 0;
for t = 2:length(xc)
   logl = logl -1/2*log(2*pi*var) - (xc(t) - xc(t-1) - theta(2)*n)^2/(2*var);
end

