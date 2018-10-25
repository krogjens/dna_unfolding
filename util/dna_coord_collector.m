function name = dna_coord_collector(expdir,pixelsize,framerate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The hp_collector routine targets a folder, specifiec by "expdir"
% and looks for trajectory files.
% The routine searches for individual dna molecule edge evolutions
% and saves these under names corresponding to the 
% experiment in the same directory as the coordinates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Target a directory with data
target = expdir;

%Specify temperature
T = 293.15; % "Room" temperature for not specified Ts

   %Find the number of independent data files (unique number before _leftHairpin for example)
   files = dir(target);
   n_traj = 0;
   data = [];
   %Find data ids
   for j = 1:length(files)
      findres =strfind(files(j).name,'leftHairpin.txt');
      isleft = max(size(findres));
      if isleft == 1 % Use the 'leftHairpin' file to extract fileid
         n_traj = n_traj + 1;
         % Add data id to data struct
         data(n_traj).fileid = files(j).name(1:(findres-1));
      end
   end
   fprintf('Found %i set(s) of trajectories in %s\n',n_traj,target);

   for j = 1:n_traj
      % Load the unfolding hairpins
      % Load left hairpin
      tarfile = [data(j).fileid,'leftHairpin.txt'];
      path = [target,'/',tarfile];
      impo_l = importdata(path);
      xl = impo_l.data(:,3) - impo_l.data(:,2);

      % Load right hairpin
      tarfile = [data(j).fileid,'rightHairpin.txt'];
      path = [target,'/',tarfile];
      impo_r = importdata(path);
      xr = impo_r.data(:,3) - impo_r.data(:,2);

      t_last = min(length(xl),length(xr));
      % Force xl to be the longer hairpin
      if length(xl) < length(xr)
         [xl,xr] = deal(xr,xl); %Swap contents of arrays
      end

      % Estimate length
      L = mean(2*(xr + xl(1:t_last)) + impo_r.data(1:t_last,2) - impo_l.data(1:t_last,3)); 

      % Save data
      data(j).xl = xl;
      data(j).xr = xr;
      data(j).L = L;
      data(j).T = T;
      data(j).deltat = 1/framerate;
      data(j).pixelsize = pixelsize;

   end
   % Name each folders data set
   name = [expdir,'/',expdir,'_data_collected'];
   save(name,'data');


