
% add all toolkit files to matlab search directory, assuming your toolkit
% is in ~/
kit_dir    = '~/DAKIT/'; % location of your analysis toolkit
path(genpath(kit_dir),path);       % add all directories in kit_dir to path

healykit_dir    = '~/Documents/MATLAB/HealyToolkit/'; % location of your analysis toolkit
path(genpath(healykit_dir),path);       % add all directories in kit_dir to path

%All the series which have the catalog waveforms
%Series which errored - Golden BH
series = cellstr(char('HR-series'))%,,'LL-series'))%,....
           %  'Lq-series', 'MaxSpin','NG-series','Oq-series','Q-series','RO1-series','RO2-series','RO3-series',....
           %  'Seff-series','SO-series','Sq-series', 'S-series-v2', 'S-series-v3', 'ssq2-opp-series', 'SSq2-series', 'SS-series', 'TP2-series', 'TP-series','A-series','EccQ','Eq-series','GoldenBH', 'GW15-series','HRq-series'))

%Path where final waveforms files will be created and stored.
numrelpath = '/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/';  %'/nethome/numrel/datafiles/Finalized_Waveforms/Waveform_files/New_BK';

%Create a temporary directory and add the path here - the strain files will
%be generated here and then will be moved to numrelpath
outdirpath = '/localdata/bkhamesra3/Test/Test'; 


% Compiling the catalog
for i=1:length(series)
     Y = sc_search('setname',series(i), 'verbose')
     for j = 76:length(Y) 

         [~,name,ext] = fileparts(Y(j).init_info.par_file_string);
         wfpath = Y(j).init_info.sim_dir;
         fprintf('%d. Creating Strain for %s waveform in %s series \n \n ',j, name,series{i})
	 try
             func_Strain(outdirpath, wfpath, numrelpath)     	
	 catch
	     j = j+1
  	 end
     end
end

%% Ignore Anything Below This

%Directories to be carefule about  - 1. MaxSpin waveforms - not all modes present and all files in .gz format. %Results with HR-series incorrect - need to sort out

%Test this later - /nethome/numrel/datafiles/Waveforms/SO-series/SO_D9_q1.5_th2_135_ph1_45_m140/SO_D9_q1.5_th2_135_ph1_45_m140-all'
%'/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.000_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.045_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.090_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.135_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.180_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.225_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.270_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.315_M77'

%wfpath = cellstr(char('/localdata2/bkhamesra3/simulations/Stampede/BBH/Event_Runs/BBH_Jan4event_UID5_M120/BBH_Jan4event_UID5_M120/D12_q2.32540_a1z_-0.12558_a2z_-0.36340_M120'))
%wfpath = cellstr(char( '/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.6/fr_b3.1_a0.6_oth.150_M77'))

%clc; clear; 
%ligoscripts_dir = '/localdata/bkhamesra3/research_localdata/UsefulScripts/LIGO/LIGO_Scripts/LIGO_Metadata/LIGO_Scripts'; 
%path(genpath(ligoscripts_dir),path);

%% Testing

% Y = sc_search('setname','D10_q5.00_a0.0_0.0_m240', 'verbose');
% wfpath = Y(1).init_info.sim_dir
% %ComputeStrain(Y(2), lmax, numrelpath)    ;	
% func_Strain( outdirpath, wfpath, numrelpath)

    
