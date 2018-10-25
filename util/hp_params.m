function [params] = hp_params(theta,M) 
n = 2; % Number of paramaters that are always active;
params=zeros(length(M)+n,1);
params(1:n)=theta(1:n);
d = 0;
if M(1) == 0
   params(3) = 1;
   d = 1;
end

for j=1+d:length(M)
  if M(j)==1 %Only transform variable parameters
    n=n+1;
    params(j+2)=theta(n);
  end
end

