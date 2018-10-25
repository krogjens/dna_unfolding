function dna_likelihood_plotter(kymo,exp_name)
l_size = 17;

C = imread([exp_name,'/',kymo(1:end-4),'.png']);

imagesc([0.2 0.6],[0.5 0.9],C)
kp = subplot('Position',[0. 0.25 0.3 .8]);
kp.XTick = [];
kp.YTick = [];
ky =subimage(C);
title('Sample kymograph','Interpreter','Latex','FontSize',l_size +2)

output = importdata([exp_name,'/',kymo(1:end-14),'_bayesout.mat']);

logl =@(obs,theta) hp_logl_lund(obs,theta,1);

model = 1;
MM = [0 0; 1 0; 2 0; 0 1; 1 1; 2 1];
for i = 1:length(output.results(model).param_mean)
    theta(i) = output.results(model).param_mean(i);
end
theta = hp_params(theta,MM(model,:));

obs = output.data;
up = 0.05;
dt = output.data.deltat;
ps = output.data.pixelsize;
temp = output.data.T;
tl = (1:length(obs.xl))*dt;
tr = (1:length(obs.xr))*dt;
tp = subplot('Position',[0.1 0.2+up 0.33 0.55]);
down = 5;
plot(tl,obs.xl*ps*10^6,'-')
hold on
plot(tr,obs.xr*ps*10^6,':','LineWidth',2)
ccs = get(gca,'colororder');
line([205 235]*.059-6,[90 90]/7,'Color',ccs(1,:))
text(240*0.059-6,90/7,'$x_l$','Interpreter','latex','FontSize',l_size)
text(240*0.059-6,80/7,'$x_r$','Interpreter','latex','FontSize',l_size)
line([205 235]*0.059-6,[80 80]/7,'LineStyle',':','LineWidth',2,'Color',ccs(2,:))

hold off
xlabel('$t \,/$ s','Interpreter','latex','FontSize',l_size)
ylabel('$x / \mu$m','Interpreter','latex','FontSize',l_size)
tp.FontSize = l_size - 4;
grid on

gamma = 0.00003*(1:170)*4;
f = 0.01*(1:100)*1.2;
mu = 0.02*(1:100);
K = 0.01*(1:100);
   gamma_phys = physconst('Boltzmann')*temp*dt/ps^3*10^3* gamma;
   f_phys = physconst('Boltzmann')*temp/ps *f *10^15;

text(15.5,22 - down -2 ,'Likelihood','Interpreter','latex','FontSize',l_size+2)
text(1,22 -down -2,'Hairpin evolution','Interpreter','latex','FontSize',l_size+2)
ll = [];
for i = 1:length(f) 
   ll(i) = logl(obs,[theta(1) f(i) theta(3) theta(4)]);
end
like = exp(ll-max(ll)); 
ax3 = subplot('Position',[0.75 0.2+up 0.20 0.55]);
plot(f_phys,like)
grid on
xlabel('$f \, /$ fN','Interpreter','latex','FontSize',l_size)
ylabel('$P(data|\theta)$','Interpreter','latex','Rotation',90,'FontSize',l_size)

ax3.YTick = [];
ll = [];
for i = 1:length(gamma) 
   ll(i) = logl(obs,[gamma(i) theta(2) theta(3) theta(4)]);
end
like = exp(ll-max(ll)); 
ax4 = subplot('Position',[0.5 0.2+up 0.20 0.55]);
plot(gamma_phys,like)
grid on
txtf = sprintf('$f = %.1f \\pm %.1f $ fN',output.results(1).param_mean(2)*physconst('Boltzmann')*temp/ps*10^15 ... 
    ,output.results(1).param_stddev(2)*(physconst('Boltzmann')*temp/ps*10^15));
txtg = sprintf('$\\gamma = %.2f \\pm %.2f $ mPa s',output.results(1).param_mean(1)*physconst('Boltzmann')*temp*dt/ps^3*10^3 ... 
    ,output.results(1).param_stddev(1)*physconst('Boltzmann')*temp*dt/ps^3*10^3);

xlabel('$\gamma \, / $ mPa $\cdot$ s','Interpreter','Latex','FontSize',l_size)
ylabel('$P(data|\theta)$','Interpreter','latex','Rotation',90,'FontSize',l_size)

text(-6,-0.3-up,txtf,'Interpreter','Latex','FontSize',16)
text(1,-0.3-up,txtg,'Interpreter','Latex','FontSize',16)
annotation('rectangle',[0.18 0.025 0.7 0.07])

ax3.FontSize=l_size-4;
ax4.FontSize=l_size-4;

hold off
fig=gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 14/2 10/2];
print([exp_name,'/',kymo(1:end-14),'_likeplot'],'-depsc')

