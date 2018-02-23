%% Create the directory of the waveform
%% Problematic waveforms - Golden Series - All file in .asc.gz form - Change the Section 2. 


wfpath = '/localdata2/bkhamesra3

lmax = 2

%if exist(dirpath, 'dir')==7
%    warning('Output directory exists at this path')
%else
%    mkdir(dirpath)
%end

[path, name, ext] = fileparts(wfpath);
%direc = fullfile(dirpath, [name, ext]);

%if exist(direc,'dir')==7
%     warning('Waveform Directory Already exist inside output directory')
%else 
%    mkdir(direc);
%end

%% Copy the necessary files from corresponding numrel - waveforms directory

shifttracker0 = fullfile(wfpath , 'ShiftTracker0.asc');
shifttracker1 = fullfile(wfpath , 'ShiftTracker1.asc');
parfile = fullfile(wfpath, [name,ext,'.par']) ;
ihspin0 = fullfile(wfpath , 'ihspin_hn_0.asc');
ihspin1 = fullfile(wfpath , 'ihspin_hn_1.asc');
ihspin3 = fullfile(wfpath , 'ihspin_hn_3.asc');
ihspin4 = fullfile(wfpath , 'ihspin_hn_4.asc');
hnmass0 = fullfile(wfpath , 'hn_mass_spin_0.asc');
hnmass1 = fullfile(wfpath , 'hn_mass_spin_1.asc');
bhdiag0 = fullfile(wfpath , 'BH_diagnostics.ah1.gp');
bhdiag1 = fullfile(wfpath , 'BH_diagnostics.ah2.gp');

%copyfile(shifttracker0, direc);
%copyfile(shifttracker1, direc);
%copyfile(parfile, direc);

%if exist(hnmass0) && exist(hnmass1)
%copyfile(hnmass0, direc);
%copyfile(hnmass1, direc);
%end 

%if exist(bhdiag0)
%    copyfile(bhdiag0, direc)
%end

%if exist(bhdiag1)
%    copyfile(bhdiag1, direc)
%end

%ihfile = 1;
%if exist(ihspin0) & exist(ihspin1)
%    copyfile(ihspin0, direc);
 %   copyfile(ihspin1, direc);
%elseif exist(ihspin3) & exist(ihspin4)
 %   copyfile(ihspin3, direc);
 %   copyfile(ihspin4, direc);
%else
 %   fprintf('Spin information not available')
 %   ihfile = 0;
%end 
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
if (ihfile==1)
   bn.init_rundata;
   sp = bn.rundata.BHp.spin;
   sm = bn.rundata.BHm.spin ;
    
end
wf = bn.calcStrain;


%% Make the Plots
n = bn.ll(end);
max_idx = n^2 + 2*n - 3;
ll = 1, mm = 1;

%plot_direc = mkdir(fullfile(direc, '/strainplots'))
% for j = 1:max_idx
%     [ll,mm] = bn.index_to_lm(j)
%     figure
%     plot(bn.time, bn.Real(:,j),'b','DisplayName', 'Real')
%     hold on
%     plot(bn.time, bn.Imag(:,j), 'g', 'DisplayName', 'Imag')
%     xlabel('Time')
%     ylabel('Strain')
%     legend('show')
%     strainplot = [plot_direc, '/', sprintf('StrainPlot_l%i_m%i.png', ll, mm)]
%     saveas(gcf, strainplot)
% end 
   

%% Output data in text file


%for idx  = 1:max_idx
%    [l,m] = bn.index_to_lm(idx);
    
%    strain_data = [direc,'/',sprintf('Strain_l%i_m%i.txt',l , m)];
%    strain_file = fopen(strain_data,'w');
%    fprintf(strain_file, '#(l,m) = %i  %i \n', l,m)
%    fprintf(strain_file, '# Time \t Real \t Imag \t Amp \t Phase \n')
%    if l==2 & m==2
%        fprintf(strain_file, '#Momega_22 = %i \n', bn.omega22start)
%        end
        
    %fprintf(strain_data, '%5d %5d %5d %5d %5d \n', wf.time(:), wf.Real(:,idx), wf.Imag(:,idx), wf.Ampl(:,idx), wf.Phse(:,idx))
%    writedata = [wf.time, wf.Real(:,idx), wf.Imag(:, idx), wf.Ampl(:,idx), wf.Phse(:,idx)];
%    dlmwrite(strain_data, writedata, 'delimiter', '\t','-append', 'precision','%0.10f');
%    fclose(strain_file);
%end

%if ihfile==1
%    spin_data = strcat(direc, '/Interpolated_Spindata.txt');
%    spin_file = fopen(spin_data,'w');
%    fprintf(spin_file, '# Time \t S1x \t S1y \t S1z \t S2x \t S2y \t S2z \n');
%    writedata = [sp.t, sp.S(:,1), sp.S(:,2), sp.S(:,3), sm.S(:,1), sm.S(:,2), sm.S(:,3)];
%    dlmwrite(spin_data, writedata, 'delimiter', '\t', 'precision','%0.10f');
%    fclose(spin_file);
%end


%if spintype == 0        %'NonSpinning'
%    numrelpath = fullfile(numrelpath, 'NonSpinning')
%elseif spintype == 1       %'AlignSpins'
%    numrelpath = fullfile(numrelpath, 'Aligned_Spins')
%elseif spintype == 2     %'Precessing'
%    numrelpath = fullfile(numrelpath, 'Precessing')
%end

%movefile(direc, numrelpath)
%clear();  
%end  
