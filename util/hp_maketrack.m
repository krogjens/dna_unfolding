function data = hp_maketrack(dat,theta,n)


gamma = theta(1);
f = theta(2);
mu = theta(3);
%var_n = theta(4);
K = theta(4);

T = 30*length(dat.xl);
L = dat.L ;
xl(1) = dat.xl(1);
xr(1) = dat.xr(1);
xb(1) = L - xl(1) - xr(1);

mean_l = -2 * f / (4 * (gamma * xl(1)^mu + K)) - ...
    mu*gamma*xl(1)^(mu-1)/4 * 1/(gamma * xl(1)^mu + K)^2; 
mean_r = -2 * f / (4 * (gamma * xr(1)^mu + K)) - ...
    mu*gamma*xr(1)^(mu-1)/4 * 1/(gamma * xr(1)^mu + K)^2; 
eta = randn(1,2);
mob(1,1) = 1/4*(1/(gamma*xl(1)^mu + K) + 1/(gamma*xb(1)^mu + 2*K));
mob(1,2) = -1/(4*(gamma*xb(1)^mu + 2*K));
mob(2,1) = mob(1,2);
mob(2,2) = 1/4*(1/(gamma*xr(1)^mu + K) + 1/(gamma*xb(1)^mu + 2*K));
[P,D] = eig(mob);

xr(2) = 0;
while ~(xr(2)> 0)
    eta =randn(1,2);
    noise = P*sqrt(2*n*D) * eta';
    xl(2) = xl(1) + mean_l*n + noise(1);
    xr(2) = xr(1) + mean_r*n + noise(2);
end
xb(2) = L - xl(2) - xr(2);
lfinish = 0;
rfinish = 0;
t = 3;

while min(rfinish,lfinish) == 0 
    mean_l = -2 * f / (4 * (gamma * xl(t-1)^mu + K)) - ...
        mu*gamma*xl(t-1)^(mu-1)/4 * 1/(gamma * xl(t-1)^mu + K)^2; 
    mean_r = -2 * f / (4 * (gamma * xr(t-1)^mu + K)) -  ...
        mu*gamma*xr(t-1)^(mu-1)/4 * 1/(gamma * xr(t-1)^mu + K)^2;
    
    eta = randn(1,2);
    mob(1,1) = 1/4*(1/(gamma*xl(t-1)^mu + K) + 1/(gamma*xb(t-1)^mu + 2*K));
    mob(1,2) = -1/(4*(gamma*xb(t-1)^mu + 2*K));
    mob(2,1) = mob(1,2);
    mob(2,2) = 1/4*(1/(gamma*xr(t-1)^mu + K) + 1/(gamma*xb(t-1)^mu + 2*K));
    [P,D] = eig(mob);
    noise = P*sqrt(2*D*n) * eta';
    xl(t) = xl(t-1) + mean_l*n + noise(1);
    xr(t) = xr(t-1) + mean_r*n + noise(2);
    xb(t) = L - xl(t) - xr(t);
   if xl(t) < 0 
       xl(t) = 0;
       lfinish = 1;
   end
   if xr(t) < 0 
       xr(t) = 0;
       rfinish = 1;
   end
   if lfinish == 1 || rfinish == 1
       tfinish = t;
       break;
   end
   t = t + 1;
end

if lfinish == 0
    obs = xl;
elseif rfinish == 0
    obs = xr;
else
    obs = zeros(1,tfinish);
end
   
if obs(tfinish) > 0 && min(rfinish,lfinish) == 0 
    t = tfinish + 1;
    while t < T
        mobb = 1/4 * (1/(gamma * obs(t-1)^mu + K) + 1/(gamma * (L - obs(t-1))^mu + K));
        var_d = 2 * mobb;
        mean_x = - 2 * f * mobb - mu*gamma / 4 *...
            ( obs(t-1)^(mu-1) / ( gamma * obs(t-1)^mu + K )^2  ...
             - (L - obs(t-1))^(mu-1) / (gamma * (L - obs(t-1))^mu + K)^2);
        var_tot = var_d * n;
        step = mean_x * n + sqrt(var_tot) * randn;
        obs(t) = obs(t-1) + step;
        if obs(t) < 0
            obs(t) = 0;
            break;
        end
        t = t + 1;
    end
end
obs = obs(~(obs==0));
if lfinish == 0
    xl = obs;
elseif rfinish == 0
    xr = obs;
end
if length(xl) < length(xr)
    data.xl = xr;
    data.xr = xl;
else
    data.xl = xl;
    data.xr = xr;
end
data.L = L;
%data.xb = xb;
%data.t_last = tfinish;

%fprintf('Generated track of %i / %i observations after %i tries\n',length(obs),T,tries);

