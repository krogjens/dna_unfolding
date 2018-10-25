function data = hp_scaling(dat,n)

n_steps = floor((length(dat.xl)-1)/n);
data.xl = dat.xl(1+(0:n_steps)*n);

if length(dat.xr) > n
    n_steps = floor((length(dat.xr)-1)/n);
    data.xr = dat.xr(1+(0:n_steps)*n);
else
    data.xr = [];
end

data.L = dat.L;
