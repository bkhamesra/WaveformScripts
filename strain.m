
% add all toolkit files to matlab search directory, assuming your toolkit
% is in ~/
kit_dir    = '~/DAKIT/'; % location of your analysis toolkit
path(genpath(kit_dir),path);       % add all directories in kit_dir to path

healykit_dir    = '~/Documents/MATLAB/HealyToolkit/'; % location of your analysis toolkit
path(genpath(healykit_dir),path);       % add all directories in kit_dir to path

%All the series which have the catalog waveforms
%Series which errored - Golden BH
series = cellstr(char('s-series'))%,,,....
           %  ,,, 'S-series-v2', 'S-series-v3', 'ssq2-opp-series', 'SSq2-series', 'SS-series', 'TP2-series', 'TP-series','A-series','EccQ','Eq-series','GoldenBH', 'GW15-series','HRq-series','HR-series', 'LL-series','Lq-series','MaxSpin','NG-series','Oq-series','Q-series','RO1-series','RO2-series','RO3-series','Seff-series','SO-series'))

%Path where final waveforms files will be created and stored.
numrelpath = '/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/';  %Path on numrel where strain output is moved

%Create a temporary directory and add the path here - the strain files will
%be generated here and then will be moved to numrelpath
outdirpath = '/localdata/bkhamesra3/Test/Test'; %Path where initial strain  output is produced


%Compiling the catalog series
% for i=1:1   %length(series)
%      Y = sc_search('setname',series(i), 'verbose')
%      for j = 4:4 %1:length(Y) 
% 
%          [~,name,ext] = fileparts(Y(j).init_info.par_file_string);
%          wfpath = Y(j).init_info.sim_dir;
%          fprintf('%d. Creating Strain for %s waveform in %s series \n \n ',j, name,series{i})
% 	     try
%              func_Strain(outdirpath, wfpath, numrelpath);
% 	     catch
% 	        j = j+1 ;
%   	     end
%      end
% end

%% Move Directories from Completed and Failed to Remaining
%remdir = '/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/AlignedSpin/';
%compdir = '/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Completed/AlignedSpin/';
%faildir = '/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Failed/AlignedSpin/';
%for i=1:length(series)
%     Y = sc_search('setname',series(i), 'verbose');
%     for j = 1:length(Y) 
%
%         [~,name,ext] = fileparts(Y(j).init_info.par_file_string);
%         wfpath = Y(j).init_info.sim_dir;
%         fprintf('%d. Creating Strain for %s waveform in %s series \n \n ',j, name,series{i});
%         
%	 wf_comppath = fullfile(compdir, name);
%	 wf_failpath = fullfile(faildir, name) ;
%	
%	 fprintf('%s \n',wf_comppath)
%	 fprintf('%s \n', wf_failpath)
%	 if(7==exist(wf_comppath, 'dir'))
%	 	movefile(wf_comppath, remdir)
%	 else
%		fprintf('Directory not found \n')
%	 end
%
%         if (7==exist(wf_failpath, 'dir')) 
%	 	movefile(wf_failpath, remdir)
%	 else
%		fprintf('Directory not found \n')
%	 end
%     end
%end
%% Ignore Anything Below This

%Directories to be carefule about  - 1. MaxSpin waveforms - not all modes present and all files in .gz format. %Results with HR-series incorrect - need to sort out

%Test this later - /nethome/numrel/datafiles/Waveforms/SO-series/SO_D9_q1.5_th2_135_ph1_45_m140/SO_D9_q1.5_th2_135_ph1_45_m140-all'
%'/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.000_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.045_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.090_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.135_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.180_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.225_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.270_M77','/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.2/fr_b3.1_a0.2_oth.315_M77'

wfpath = cellstr(char('/localdata2/bkhamesra3/simulations/Stampede/BBH/Event_Runs/BBH_Jan4event_UID5_M120/BBH_Jan4event_UID5_M120/D12_q2.32540_a1z_-0.12558_a2z_-0.36340_M120'))
%wfpath = cellstr(char( '/nethome/numrel/datafiles/Waveforms/archive/s-series/spin0.6/fr_b3.1_a0.6_oth.150_M77'))

%clc; clear; 
%ligoscripts_dir = '/localdata/bkhamesra3/research_localdata/UsefulScripts/LIGO/LIGO_Scripts/LIGO_Metadata/LIGO_Scripts';%Add path of this file  
%path(genpath(ligoscripts_dir),path);
%
%%% Testing
%
% Y = sc_search('setname','D10_q5.00_a0.0_0.0_m240', 'verbose');%Change the simulation name as desired
% wfpath = Y(1).init_info.sim_dir
wfpath =   '/nethome/numrel/datafiles/Waveforms/NG-series/D9_q5.0_a0.0_Q24';
 func_Strain( outdirpath, wfpath, numrelpath)
%
    
