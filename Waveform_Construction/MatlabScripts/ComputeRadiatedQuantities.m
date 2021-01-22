%% Script to compute radiated energy, momenta and angular momenta from BBH/BHNS simulation
warning('off','all')

% Add path of HealyToolkit
healykit_dir    = '~/Documents/MATLAB/HealyToolkit/'; % location of your analysis toolkit
path(genpath(healykit_dir),path);       % add all directories in kit_dir to path

tic;

% simulation directory
simdir   = '/localdata2/bkhamesra3/simulations/Hive/NSBH/MassRatio_5/CompactnessStudy/NSBH_q5_D9_a0_C0.2_M74/Summary/Data'; 


% output directory - Path to save the radiated quantities data
strdir = fullfile(simdir, 'Strain_FiniteRadius');
if ne(exist(strdir, 'dir'),7)
    mkdir(strdir);
end

% Compute system characteristics from simulation name - Assumes form <system_type>_q<mass_ratio>_D<separarion>_a<spinmag>_C<compactness>_M<resolution
simdir_parts = strsplit(simdir, '/');
simname = char(simdir_parts(length(simdir_parts)-2));
simname_parts = strsplit(simname,'_');
system = char(simname_parts(1));
mass_ratio = char(simname_parts(2));
mass_ratio = str2double(mass_ratio(2));
comp = char(simname_parts(5));
comp = str2double(comp(2:end))

% Select the extraction radius
r_ex = [70,100]; % [50,60,70,80,90,100,110,120,130,140];

% Add the adm mass
admmass_nsbh = struct('Q5_C012',0.994416, 'Q5_C015',0.994092, 'Q5_C018',0.995021, 'Q5_C020',0.994795, 'Q3_C012',0.992497, 'Q3_C015',0.992126, 'Q3_C018',0.992665, 'Q3_C020',0.993399);

if mass_ratio==2 & strcmp(system,'NSBH') & comp==0.15
    admmass_init = 0.9913630194925424632 ;
elseif mass_ratio==3 & strcmp(system,'NSBH') & comp==0.12
    admmass_init = admmass_nsbh.Q3_C012;
elseif mass_ratio==3 & strcmp(system,'NSBH') & comp==0.15
    admmass_init = admmass_nsbh.Q3_C015;
elseif mass_ratio==3 & strcmp(system,'NSBH') & comp==0.175
    admmass_init = admmass_nsbh.Q3_C018;
elseif mass_ratio==3 & strcmp(system,'NSBH') & comp==0.2
    admmass_init = admmass_nsbh.Q3_C020;
elseif mass_ratio==5 & strcmp(system,'NSBH') & comp==0.12
    admmass_init =  admmass_nsbh.Q5_C012; 
elseif mass_ratio==5 & strcmp(system,'NSBH') & comp==0.15
    admmass_init =  admmass_nsbh.Q5_C015; %0.9940921134049298669;
elseif mass_ratio==5 & strcmp(system,'NSBH') & comp==0.175
    admmass_init =  admmass_nsbh.Q5_C018; %0.9940921134049298669;
elseif mass_ratio==5 & strcmp(system,'NSBH') & comp==0.2
    admmass_init =  admmass_nsbh.Q5_C020; %0.9940921134049298669;
elseif mass_ratio==2 & strcmp(system,'BBH') 
    admmass_init = 0.990134716347218502 ;
elseif mass_ratio==3 & strcmp(system,'BBH') 
    admmass_init = 0.9917528543336998625 ;
elseif mass_ratio==5 & strcmp(system,'BBH') 
    admmass_init = 0.9939733957185699076 ;
end

for idx = 1:length(r_ex)
    % create modeBundle object - Use FFI method 
    % Caution -  Do not extrapolate to infinity as it leads to difference in the final values of Angular momenta 
    wf_l2 = modeBundle(simdir,2:2,r_ex(idx),'Madm',admmass_init)
    wf_l3 = modeBundle(simdir,2:3,r_ex(idx),'Madm',admmass_init)
    wf_l4 = modeBundle(simdir,2:4,r_ex(idx),'Madm',admmass_init)
    wf_l5 = modeBundle(simdir,2:5,r_ex(idx),'Madm',admmass_init)
    wf_l6 = modeBundle(simdir,2:6,r_ex(idx),'Madm',admmass_init)
    
    wf    = modeBundle(simdir, 2:8, r_ex(idx), 'Madm', admmass_init)

    % Compute Initial frequency
    wf.calcOmega

    % Compute energy, momenta and angular momenta using radAway class
    rA = radAway(wf, 'rmics',1);

    %output the data - in simulation directory
    rA.output_timeseries(simdir);    

    %Save strain and its time derivative
    hdot =  rA.Hdot; 
    strain = rA.Strn; %wf.calcStrain

    max_idx = n^2 + 2*n - 3;
    for lm_idx  = 1:max_idx
        [l,m] = wf.index_to_lm(lm_idx);
        
        hdot_data = [strdir,'/',sprintf('Hdot_l%i_m%i_r%i.txt', l, m, r_ex(idx))];
    
        hdot_file = fopen(hdot_data,'w');
        fprintf(hdot_file, '#(l,m,r_ex) = %i  %i %f\n', l,m,r_ex(idx))
        fprintf(hdot_file, '# Time \t Real \t Imag \t Amp \t Phase \n')
        writedata = [hdot.time, hdot.Real(:,lm_idx), hdot.Imag(:, lm_idx), hdot.Ampl(:,lm_idx), hdot.Phse(:,lm_idx)];
        dlmwrite(hdot_data, writedata, 'delimiter', '\t','-append', 'precision','%.10g');
        fclose(hdot_file);
            
            
        strain_data = [strdir,'/',sprintf('Strain_l%i_m%i_r%i.txt', l, m, r_ex(idx))];
    
        strain_file = fopen(strain_data,'w');
        fprintf(strain_file, '#(l,m,r_ex) = %i  %i %f\n', l,m,r_ex(idx))
        fprintf(strain_file, '# Time \t Real \t Imag \t Amp \t Phase \n')
        writedata = [strain.time, strain.Real(:,lm_idx), strain.Imag(:, lm_idx), strain.Ampl(:,lm_idx), strain.Phse(:,lm_idx)];
        dlmwrite(strain_data, writedata, 'delimiter', '\t','-append', 'precision','%.10g');
        fclose(strain_file);
            
    end

    % Compute radiated quantities for various l modes and save the data (rmics=1 uses FFI method)
    rA_l2 = radAway(wf_l2, 'rmics',1);
    rA_l2.output_timeseries(simdir);    

    rA_l3 = radAway(wf_l3, 'rmics',1);
    rA_l3.output_timeseries(simdir);    

    rA_l4 = radAway(wf_l4, 'rmics',1);
    rA_l4.output_timeseries(simdir);    

    rA_l5 = radAway(wf_l5, 'rmics',1);
    rA_l5.output_timeseries(simdir);    

    rA_l6 = radAway(wf_l6, 'rmics',1);
    rA_l6.output_timeseries(simdir);    

    % calculate skymap
    tindex = 1000;
    [TH,PH] = skyMath.setup_mesh_grid( pi/45, pi/45 );
    Edot = rA.calcEPdotOverSky( tindex, TH, PH, 'E' );
    
    %% plot skymap
    plottitle=sprintf('t = %g M', rA.Psi4.time(tindex) );
    fig = figure(1); clf;
    set(fig, 'Visible', 'off')
    Edot.plot('colormap');
    colorbar;
    title(plottitle);
    %saveas(fig, fullfile(figdir,sprintf('Edot_NSBH_r%d.png', r_nsbh(idx)))); 
    close(fig);
    
    
    % plot totals
    rA_l6.print_totals();

end
toc;

