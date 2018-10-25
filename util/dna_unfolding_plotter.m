function dna_unfolding_plotter(kymo,kymodir)

clf
% Analyze all files

name = kymo(1:end-4);
path = kymodir;
%import image
C = imread([name,'.tif']);
cpath = path;
figure(1)
plot(1,1)
imshow(C,'InitialMagnification','fit');
hold on

%import trajectories
dsfile = importdata([cpath,'/',name,'_intactDsDNA.txt']);  
spRowsDsDNALeft = dsfile.data(:,1);
spColsDsDNALeft = dsfile.data(:,2);
spRowsDsDNARight = dsfile.data(:,1);
spColsDsDNARight = dsfile.data(:,3);

lhpfile = importdata([cpath,'/',name,'_leftHairpin.txt']);
t0 = lhpfile.data(1,1);
spColsLeftHairpinLeftEdge = lhpfile.data(:,2);
spRowsLeftHairpinLeftEdge = lhpfile.data(:,1);
spColsLeftHairpinRightEdge = lhpfile.data(:,3);
spRowsLeftHairpinRightEdge = lhpfile.data(:,1);
rhpfile = importdata([cpath,'/',name,'_rightHairpin.txt']);
spColsRightHairpinLeftEdge = rhpfile.data(:,2);
spRowsRightHairpinLeftEdge = rhpfile.data(:,1);
spColsRightHairpinRightEdge = rhpfile.data(:,3);
spRowsRightHairpinRightEdge = rhpfile.data(:,1);
dd = 0.05;

ccs = get(gca,'colororder');
ccs = ccs+dd;
% 1. Left hairpin, left edge trajectory
plot(spColsLeftHairpinLeftEdge,spRowsLeftHairpinLeftEdge,'Color',ccs(2,:),'LineWidth',3)
% 2. Left hairpin, right edge trajectory
plot(spColsLeftHairpinRightEdge,spRowsLeftHairpinRightEdge,'Color',ccs(2,:),'LineWidth',3)
% 3. Right hairpin, left edge trajectory
plot(spColsRightHairpinLeftEdge,spRowsRightHairpinLeftEdge,'Color',ccs(1,:),'LineWidth',3)
% 4. Right hairpin, left edge trajectory
plot(spColsRightHairpinRightEdge,spRowsRightHairpinRightEdge,'Color',ccs(1,:),'LineWidth',3)

% Mark cut time with white dashes
line([spColsLeftHairpinLeftEdge(1)-20 spColsRightHairpinRightEdge(1)+20],[t0 t0],'Color',[0.99 1 1],'LineStyle','--','LineWidth',.8)
text(spColsLeftHairpinLeftEdge(1)-50, t0,'$t = 0$','Color',[0.99 1 1],'Interpreter','Latex','FontSize',17)
dt = 30;

ylabel('Time','FontSize',17,'Interpreter','Latex')
annotation('arrow',[.105 .105],[0.7 0.35])
ylim([150 270])
xlim([70 300])
title('Kymograph and detected edges','Interpreter','Latex','FontSize',19);
hold off

fig=gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 14/2 7/2];

print([path,'/',name],'-dpng')


