function [theta] = hp_invprior(u,ranges,M)
%jeff=@(u,thmin,thmax) thmin*exp(u*log(thmax/thmin));
%uni=@(u,theta_min,theta_max) (theta_max - theta_min)*u +theta_min;
theta=zeros(length(u),1);
theta(1)=util_jeff(u(1),ranges(1,1),ranges(1,2));
theta(2)=util_uni(u(2),ranges(2,1),ranges(2,2));
n=2; % Number of parameters that are always active

for j=1:length(M)
  if M(j)==1 %Only transform variable parameters
    n = n + 1; %Update number of active parameters
    if ranges(j+2,1) > 0
       theta(n)=util_jeff(u(n),ranges(j+2,1),ranges(j+2,2));
    else
       theta(n)=util_uni(u(n),ranges(j+2,1),ranges(j+2,2));
    end
  end
end

