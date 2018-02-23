
%% Create the outdirpathtory of the waveform

function [] = ComputeStrain(Y, lmax, numrelpath)

%% Original

spin_plus = Y.init_info.s_plus
spin_minus = Y.init_info.s_minus

if any(spin_plus(1:2)) || any(spin_minus(1:2))
    spintype = 2;       %'Precessing';
    numrelpath = fullfile(numrelpath, 'Precessing');
    fprintf('This is a precessing BBH simulation \n')
elseif any(spin_plus) || any(spin_minus)
    spintype = 1;       %'AlignSpins';
    numrelpath = fullfile(numrelpath, 'AlignedSpin');
    fprintf('This is Aligned-Spin BBH simulation \n')
else
    spintype = 0 ;      %'NonSpinning'
    numrelpath = fullfile(numrelpath, 'NonSpinning');
    fprintf('This is Non-Spinning BBH simulation \n')
    
end
wfdirpath = Y.init_info.sim_dir;
outdirpath = numrelpath;

% This is temporary fix, Fix DAKIT init_info.sim_dir
if exist(wfdirpath)
    fprintf('wfdirpath exists')
else
    wfdirpath = strcat('/',Y.init_info.sim_dir)
end


[~, name, ext] = fileparts(Y.init_info.par_file_string);
outdirpath = fullfile(outdirpath, name);
fprintf('Path of Waveform = %s, and Path of output directory = %s  \n',wfdirpath, outdirpath)

if exist(outdirpath,'dir')==7
    warning('Waveform output directory already exists ')
else 
    mkdir(outdirpath);
end


%% Copy the necessary files from corresponding numrel - waveforms outdirpathtory

shifttracker0 = fullfile(wfdirpath , 'ShiftTracker0.asc');
shifttracker1 = fullfile(wfdirpath , 'ShiftTracker1.asc');
parfile = fullfile(wfdirpath, Y.init_info.par_file_string) ;
ihspin0 = fullfile(wfdirpath , 'ihspin_hn_0.asc');
ihspin1 = fullfile(wfdirpath , 'ihspin_hn_1.asc');
ihspin3 = fullfile(wfdirpath , 'ihspin_hn_3.asc');
ihspin4 = fullfile(wfdirpath , 'ihspin_hn_4.asc');
hnmass0 = fullfile(wfdirpath , 'hn_mass_spin_0.asc');
hnmass1 = fullfile(wfdirpath , 'hn_mass_spin_1.asc');
bhdiag0 = fullfile(wfdirpath , 'BH_diagnostics.ah1.gp');
bhdiag1 = fullfile(wfdirpath , 'BH_diagnostics.ah2.gp');

copyfile(shifttracker0, outdirpath);
copyfile(shifttracker1, outdirpath);
copyfile(parfile, outdirpath);
 
 if exist(hnmass0) && exist(hnmass1)
    copyfile(hnmass0, outdirpath);
    copyfile(hnmass1, outdirpath);
 end 
 
 if exist(bhdiag0)
    copyfile(bhdiag0, outdirpath)
 end
 
 if exist(bhdiag1)
    copyfile(bhdiag1, outdirpath)
 end
 
 ihfile = 1;
 if exist(ihspin0) & exist(ihspin1)
     copyfile(ihspin0, outdirpath);
     copyfile(ihspin1, outdirpath);
 end 
 
 if exist(ihspin3) & exist(ihspin4)
     copyfile(ihspin3, outdirpath);
     copyfile(ihspin4, outdirpath);
 elseif ~exist(ihspin0) & ~exist(ihspin3)
     fprintf('Spin information not available')
     ihfile = 0;
     end 
 %% Compute Strain - 
 
 figdirpath = fullfile(outdirpath,'figures')
 fprintf('figures directory - %s',figdirpath)
 
 if exist(figdirpath,'dir')==7
    warning('Waveform figures directory already exists at path - %s\n',figdirpath)
 else 
    mkdir(figdirpath);
    fprintf('Figures directory created')
 end
 
 for l = 2:lmax
     for m = -l:l
         y = y_load(Y,'lm',[l m],'rscale','verbose');
         h = y_strain(y,'safe_opt','shift',0); % note that this shift is the same as above -- y_strain calls y_clean, and the functions pass input argments to each other
         fig = figure; box on;hold on;
         plot(h.t_raw,h.Yp_raw,'Color','b')
         plot(h.t_raw,h.Yx_raw,'--r')
            saveas(fig, fullfile(figdirpath,sprintf('Strain_%d%d.png',l,m)));
         hold off;
     end
 end


end


