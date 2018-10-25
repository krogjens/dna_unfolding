function dna_bayes_skel(exp_name,inputname)

input = importdata([inputname,'.mat']);
misc.data_id = exp_name;
fprintf('Found %i data series in the %s set\n',length(input),inputname)
for j = 1:length(input)
   datid = input(j).fileid;
   fprintf('Analysing %s data within the %s set\n',datid,inputname)
   %load data 
   
   data = struct;
   data.xl = input(j).xl;
   data.xr = input(j).xr;
   data.L = input(j).L;
   data.T = input(j).T;
   data.pixelsize = input(j).pixelsize;
   data.deltat = input(j).deltat;

   %Specify prior ranges
   gammamin=10^(-5);
   gammamax=10;   

   fmin = 0;
   fmax = 20;

   mumin = 0;
   mumax = 2;

   Kmin=0; 
   Kmax=2;

   %Specify prior ranges
   ranges=[gammamin gammamax ; fmin fmax; mumin mumax; Kmin Kmax];

   %Tmin = 2; % Required minimum length of a scaled data set 
   nmax = 3; %floor((length(data.xr)-1) / Tmin)
   if ~(nmax > 1)
      fprintf('Data set too short to be scaled to length >= %i',Tmin);
   elseif nmax > 10
      nmax = 10;
   end

   %Specify options
   options.nwalkers=300;   % Number of walkers to be generated
   options.stoprat=10^(-5);% Ratio to stop sampling
   options.nsteps=40;      % Attempted number of steps in parameter space
   options.nlist = 1:1:nmax;    % Maximum time scaling factor for data
   options.trackmax = 100; % Number of replicated tracks to compare information for

   %Specify the models
%   MM=[0 0; 1 0; 2 0; 0 1; 1 1; 2 1];
	MM = [0 0]; % This model has friction \propto x^-1 like the one described in the paper.
	% In general we set friction = 1/(x^mu + K) where MM(1) = 0,2,1 makes mu = 1,0 and nonfixed respectively;
	% If M(2) == 0, then K = 0 and M(2) == 1 allows for nonzero (positive) values of K;
   n_perm = 2; % Number of permanently active parameters
   for i=1:size(MM,1)
     models(i).options=options;
     models(i).genu=@() rand(1,sum(MM(i,:) == 1)+n_perm);
     models(i).logl=@(obs,theta) hp_logl_lund(obs,hp_params(theta,MM(i,:)),1);
     models(i).logl_n=@(obs,theta,n) hp_logl_lund(obs,hp_params(theta,MM(i,:)),n);
     models(i).scaling =@(obs,n) hp_scaling(obs,n); % Function that scales data
     models(i).replicate =@(obs,theta,n) hp_maketrack(obs,hp_params(theta,MM(i,:)),n);
     models(i).invprior=@(u) hp_invprior(u,ranges,MM(i,:));
     models(i).labels = 1:n_perm;
     for j=1:length(MM(i,:))
       if MM(i,j)==1
         models(i).labels=[models(i).labels j+n_perm];
       end
     end
     models(i).add{1}=@(theta) theta(1)*physconst('Boltzmann')*data.T*data.deltat ...
         / (data.pixelsize)^3;
     models(i).add{2}=@(theta) theta(2)*physconst('Boltzmann')*data.T ...
         / (data.pixelsize);
     models(i).add{3}=@(theta) theta(2)/theta(1) * data.pixelsize^2/data.deltat;
     models(i).labels=[models(i).labels 5 6 7];
   end   
   
   %Percentiles to be calculated
   misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];
   
   %Labels for the parameters
   misc.labels=...
   ['gamma:          ';...
    'force:          ';...
    'exponent:       ';...
    'perm. friction: ';...
    'SI gamma        ';...
    'SI force        ';...
    'SI force/gamma  '];
   
   misc.titles=...
   {'Percentile at:  '};

   %Tell ns_print to write a summary-file
   misc.nssummary=['/',exp_name,'_results.txt'];
   
   %Run the nested sampling algorithm for all models and compile results
   [results] = ns_processdataset(data,models,misc);
     
   fprintf('Saving results\n')
   path=[exp_name,'/',exp_name,'_bayesout'];

   %Save results
   save(path,'results')
   save(path,'data','-append')
   save(path,'options','-append')
   save(path,'ranges','-append')
end

