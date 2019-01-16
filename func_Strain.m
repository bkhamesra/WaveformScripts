%% Create the directory of the waveform
%% Problematic waveforms - Golden Series - All file in .asc.gz form - Change the Section 2. 


function [] = func_Strain(outdirpath, wfpath, numrelpath)


% add all toolkit files to matlab search directory, assuming your toolkit
% is in ~/
if exist(outdirpath, 'dir')==7
    warning('Output directory exists at this path')
else
    mkdir(outdirpath)
end

%[path, name, ext] = fileparts(wfpath);
%direc = fullfile(dirpath, [name, ext]);
wfpath_parts = strsplit(wfpath, '/');
if strcmp(wfpath_parts(length(wfpath_parts)),'')==1
    wfpath_parts = wfpath_parts(1:length(wfpath_parts)-1);
    wfpath = strjoin(wfpath_parts,'/')    ;

end    
   
wfpath_parts = strsplit(wfpath, 'nethome');

if strcmp(wfpath_parts(1),'/')==0
    wfpath = strcat('/',wfpath);
end

dirname_parts = strsplit(wfpath, '/');
dirname = dirname_parts(length(dirname_parts));
outdir =fullfile(outdirpath, strcat(dirname{1}, '/'));


if exist(outdir,'dir')==7
     warning('Waveform Directory Already exist inside output directory')
else 
     mkdir(outdir);
end

datadir = fullfile(outdir, 'data');
figdir = fullfile(outdir, 'figures');
strdir = fullfile(datadir, 'Strain');

fprintf('%f \n', exist(datadir,'dir'))
if (exist(datadir,'dir')==0)
    mkdir(datadir)
end

if ne(exist(figdir,'dir'),7)
    mkdir(figdir)
end
   
if ne(exist(strdir, 'dir'),7)
    mkdir(strdir)
end


%% Copy the necessary files from corresponding numrel - waveforms directory

%dirname_check = strsplit(dirname{1}, '-');
%if strcmp(dirname_check(length(dirname_check)),'all')
%    simname = strcat(dirname_check(1:length(dirname_check)-1));
%else
%    simname = dirname{1};
%end

simname = strsplit(dirname{1},'-all')
fprintf('%s \n',simname{1})
parfile = fullfile(wfpath,strcat(simname{1}, '.par'))    % [dirname{1},'.par']) ;
shifttracker0 = fullfile(wfpath , 'ShiftTracker0.asc');
shifttracker1 = fullfile(wfpath , 'ShiftTracker1.asc');
ihspin0 = fullfile(wfpath , 'ihspin_hn_0.asc');
ihspin1 = fullfile(wfpath , 'ihspin_hn_1.asc');
ihspin3 = fullfile(wfpath , 'ihspin_hn_3.asc');
ihspin4 = fullfile(wfpath , 'ihspin_hn_4.asc');
hnmass0 = fullfile(wfpath , 'hn_mass_spin_0.asc');
hnmass1 = fullfile(wfpath , 'hn_mass_spin_1.asc');
bhdiag0 = fullfile(wfpath , 'BH_diagnostics.ah1.gp');
bhdiag1 = fullfile(wfpath , 'BH_diagnostics.ah2.gp');
psi4_l8 = fullfile(wfpath , 'Ylm_WEYLSCAL4::Psi4r_l8_m8_r75.00.asc'); 
psi4_l6 = fullfile(wfpath , 'Ylm_WEYLSCAL4::Psi4r_l6_m6_r75.00.asc') ;
mp_l8 = fullfile(wfpath, 'mp_WeylScal4::Psi4i_l8_m8_r75.00.asc');
mp_l6 = fullfile(wfpath, 'mp_WeylScal4::Psi4i_l6_m6_r75.00.asc');

copyfile(shifttracker0, datadir);
copyfile(shifttracker1, datadir);
copyfile(parfile, datadir)

if exist(hnmass0) && exist(hnmass1)
copyfile(hnmass0, datadir);
copyfile(hnmass1, datadir);
end 

if exist(bhdiag0)
    copyfile(bhdiag0, datadir);
end

if exist(bhdiag1)
    copyfile(bhdiag1, datadir);
end

ihfile = 1;
if exist(ihspin0) & exist(ihspin1)
    copyfile(ihspin0, datadir);
    copyfile(ihspin1, datadir);
elseif exist(ihspin3) & exist(ihspin4)
    copyfile(ihspin3, datadir);
    copyfile(ihspin4, datadir);
else
    fprintf('Spin information not available')
    ihfile = 0;
end 

%% Find lmax
if (exist(psi4_l8,'file')==2)
    lmax = 8;
elseif (exist(mp_l8, 'file')==2)
    lmax = 8;
elseif (exist(psi4_l6,'file')==2) 
    lmax = 6;
elseif (exist(mp_l6, 'file')==2)
    lmax = 6;
else 
    lmax = 2;
end

lmax=2
fprintf('lmax =  %d \n', lmax)
%% Check precessing/non-precessing
spin_plus = [];

for k = 0:2
    search_txt = sprintf('"twopunctures::par_s_plus%s%i%s*"','\[',k,'\]');
    command = sprintf('grep -r %s %s',search_txt , parfile);
    [status,cmdout] = system(command);
    sp = strsplit(cmdout, '=');
    sp = str2double(strtrim(sp(length(sp))));
    spin_plus = [spin_plus,sp];
end

spin_minus = [];

for k = 0:2
    search_txt = sprintf('"twopunctures::par_s_minus%s%i%s*"','\[',k,'\]');
    command = sprintf('grep -r %s %s',search_txt , parfile);
    [status,cmdout] = system(command);
    sm = strsplit(cmdout, '=');
    sm = str2double(strtrim(sm(length(sm))));
    spin_minus = [spin_minus,sm];
end

if any(spin_plus(1:2)) || any(spin_minus(1:2))
    spintype = 2;       %'Precessing';
    fprintf('This is a precessing BBH simulation \n')
elseif any(spin_plus) || any(spin_minus)
    spintype = 1;       %'AlignSpins';
    fprintf('This is Aligned-Spin BBH simulation \n')
else
    spintype = 0 ;      %'NonSpinning'
    fprintf('This is Non-Spinning BBH simulation \n')
    
end
      


%% Create the strain data

l = 2:lmax  ;               %Change if all modes not present
r = 75;

bn = modeBundle(wfpath, l, r)
bn.calcOmega
bn.extrapolate(2)
if (ihfile==1)
   bn.init_rundata;
   sp = bn.rundata.BHp.spin;
   sm = bn.rundata.BHm.spin ;
    
end
wf = bn.calcStrain;
wf.window_ends();



        
%% Make the Plots
n = bn.ll(end);
max_idx = n^2 + 2*n - 3;
ll = 1, mm = 1;


for j = 1:max_idx
    [ll,mm] = bn.index_to_lm(j);
    fig = figure;
    set(fig, 'Visible', 'off')
    plot(wf.time, wf.Real(:,j),'b','DisplayName', 'Real');
    hold on;
    plot(wf.time, wf.Imag(:,j), 'k', 'DisplayName', 'Imag');
    xlabel('Time');
    ylabel('Strain');
    legend('show');
    saveas(fig, fullfile(figdir,sprintf('Strain_%d%d.png',ll,mm)));
    
end 
   

if (spintype==2)
    for j = 1:max_idx
        [ll,mm] = bn.index_to_lm(j);
        if ((ll==2) && (mm==2))
            mm_plus2 =  j;
        elseif((ll==2) && (mm==-2))
            mm_minus2 = j;
            
        end
    end
    
    fprintf('m=2 index = %d, m=-2 index = %d \n',mm_plus2,mm_minus2); 
    fig = figure;
    set(fig, 'Visible', 'off')
    plot(wf.time, wf.Real(:,mm_plus2),'k', 'DisplayName', 'Real');
    hold on;
    plot(wf.time, wf.Real(:,mm_minus2), '--b', 'DisplayName', 'Real');
    hold off;
    xlabel('Time');
    ylabel('Strain Real');
    legend('show');
    saveas(fig, fullfile(figdir,sprintf('Strain_l2m2_l2m-2_realcomp.png')));
    
    fig = figure;
    set(fig, 'Visible', 'off')
    plot(wf.time, wf.Imag(:,mm_plus2),'k', 'DisplayName', 'm=2');
    hold on;
    plot(wf.time, wf.Imag(:,mm_minus2), '--b', 'DisplayName', 'm=-2');
    hold off;
    xlabel('Time');
    ylabel('Strain Imag');
    legend('show');
    saveas(fig, fullfile(figdir,sprintf('Strain_l2m2_l2m-2_imagcomp.png')));
    
    fig = figure;
    set(fig, 'Visible', 'off')
    plot(wf.time, wf.Ampl(:,mm_plus2),'k', 'DisplayName', 'm=2');
    hold on;
    plot(wf.time, wf.Ampl(:,mm_minus2), '--b', 'DisplayName', 'm=-2');
    hold off;
    xlabel('Time');
    ylabel('Strain Imag');
    legend('show');
    saveas(fig, fullfile(figdir,sprintf('Strain_l2m2_l2m-2_ampcomp.png')));
    
end
%% Output data in text file


for idx  = 1:max_idx
    [l,m] = bn.index_to_lm(idx);
    
    strain_data = [strdir,'/',sprintf('Strain_l%i_m%i.txt',l , m)];

    strain_file = fopen(strain_data,'w');
    fprintf(strain_file, '#(l,m) = %i  %i \n', l,m)
    fprintf(strain_file, '# Time \t Real \t Imag \t Amp \t Phase \n')
    if l==2 & m==2
        fprintf(strain_file, '#Momega_22 = %i \n', bn.omega22start)
        end
        
    %fprintf(strain_data, '%5d %5d %5d %5d %5d \n', wf.time(:), wf.Real(:,idx), wf.Imag(:,idx), wf.Ampl(:,idx), wf.Phse(:,idx))
    writedata = [wf.time, wf.Real(:,idx), wf.Imag(:, idx), wf.Ampl(:,idx), wf.Phse(:,idx)];
    dlmwrite(strain_data, writedata, 'delimiter', '\t','-append', 'precision','%0.10f');
    fclose(strain_file);
end

%Spin Interpolated Data - Does not give correct results for various
%waveforms
% if ihfile==1
%     spin_data = strcat(direc, '/Interpolated_Spindata.txt');
%     spin_file = fopen(spin_data,'w');
%     fprintf(spin_file, '# Time \t S1x \t S1y \t S1z \t S2x \t S2y \t S2z \n');
%     writedata = [sp.t, sp.S(:,1), sp.S(:,2), sp.S(:,3), sm.S(:,1), sm.S(:,2), sm.S(:,3)];
%     dlmwrite(spin_data, writedata, 'delimiter', '\t', 'precision','%0.10f');
%     fclose(spin_file);
% end


if spintype == 0        %'NonSpinning'
    numrelpath = fullfile(numrelpath, 'NonSpinning')
elseif spintype == 1       %'AlignSpins'
    numrelpath = fullfile(numrelpath, 'AlignedSpin')
elseif spintype == 2     %'Precessing'
    numrelpath = fullfile(numrelpath, 'Precessing')
end

movefile(outdir, numrelpath)
clear();  
end  
