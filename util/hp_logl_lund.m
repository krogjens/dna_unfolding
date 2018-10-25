function logl = hp_logl_lund(data,theta,n)

% obs must be a Tx3 array with coordinates for the two hairpin lengths 
% and the DNA center, respectively. xl is the hairpin that survives the
% longest;

gamma = theta(1);
f = theta(2);
mu = theta(3);
%var_n = theta(4);
K = theta(4);
xl = data.xl;
xr = data.xr;
L = data.L; 
xb = L - xl(1:length(xr)) - xr;
%xb = data.xb;
t_last = max(length(xr) - 1,0);
T = length(xl);

logl = 0;
if t_last > 0
    mean_l = -2 * f / (4 * (gamma * xl(1)^mu + K)) - ...
        mu*gamma*xl(1)^(mu-1)/4 * 1/(gamma * xl(1)^mu + K)^2;
    mean_r = -2 * f / (4 * (gamma * xr(1)^mu + K)) - ...
        mu*gamma*xr(1)^(mu-1)/4 * 1/(gamma * xr(1)^mu + K)^2;
    y=zeros(1,2);
    y(1) = xl(2) - xl(1) - mean_l*n;
    y(2) = xr(2) - xr(1) - mean_r*n;
    det = (4*K + gamma*(xb(1)^mu + xl(1)^mu + xr(1)^mu)) ...
        /(16*(gamma * xl(1)^mu + K)*(gamma * xr(1)^mu + K)*(gamma * xb(1)^mu + 2*K));
    
    logl = logl - log(4*pi*n) - 1/2 *log(det) - 1/(4 * n) * 1/(4*det) * ...
        (y(1)^2/(gamma*xr(1)^mu + K) + y(2)^2/(gamma*xl(1)^mu + K) + (y(1) + y(2))^2/(gamma*xb(1)^mu + 2*K));
    
    for t = 3:t_last + 1
        mean_l = -2 * f / (4 * (gamma * xl(t-1)^mu + K)) - ...
            mu*gamma*xl(t-1)^(mu-1)/4 * 1/(gamma * xl(t-1)^mu + K)^2;
        mean_r = -2 * f / (4 * (gamma * xr(t-1)^mu + K)) - ...
            mu*gamma*xr(t-1)^(mu-1)/4 * 1/(gamma * xr(t-1)^mu + K)^2;
        y=zeros(1,2);
        y(1) = xl(t) - xl(t-1) - mean_l*n;
        y(2) = xr(t) - xr(t-1) - mean_r*n;
        det = (4*K + gamma*(xb(t-1)^mu + xl(t-1)^mu + xr(t-1)^mu)) ...
            /(16*(gamma * xl(t-1)^mu + K)*(gamma * xr(t-1)^mu + K)*(gamma * xb(t-1)^mu + 2*K));
        logl = logl - log(4*pi*n) - 1/2 *log(det) - 1/(4 * n) * 1/(4*det) * ...
            (y(1)^2/(gamma*xr(t-1)^mu + K) + y(2)^2/(gamma*xl(t-1)^mu + K) + (y(1) + y(2))^2/(gamma*xb(t-1)^mu + 2*K));
    end
end

for t = t_last + 2:T
   mob = 1/4 * (1/(gamma * xl(t-1)^mu + K) + 1/(gamma * (L - xl(t-1))^mu + K)); 
   var_d = 2 * mob * n;
   mean_x = - 2 * f * mob - mu *gamma/4 * ...
      (xl(t-1)^(mu-1)/(gamma*xl(t-1)^mu + K)^2 ... 
        -(L-xl(t-1))^(mu-1)/(gamma*(L-xl(t-1))^mu + K)^2);
   F = xl(t) - xl(t-1) - mean_x*n;
   var_tot = var_d;
   logl = logl - log(sqrt(2*pi*var_tot)) - F^2 /(2 * var_tot);
end

    
