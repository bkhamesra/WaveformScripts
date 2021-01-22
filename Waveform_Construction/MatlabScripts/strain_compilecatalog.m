clc; clear;

% add all toolkit files to matlab search directory, assuming your toolkit
% is in ~/
kit_dir    = '~/DAKIT/'; % location of your analysis toolkit
path(genpath(kit_dir),path);       % add all directories in kit_dir to path

healykit_dir    = '~/Documents/MATLAB/HealyToolkit/'; % location of your analysis toolkit
path(genpath(healykit_dir),path);       % add all directories in kit_dir to path


%Path where final waveforms files will be created and stored.
numrelpath = '/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/';  %Path on numrel where strain output is moved

%Create a temporary directory and add the path here - the strain files will
%be generated here and then will be moved to numrelpath
outdirpath = '/localdata2/bkhamesra3/Waveforms'; %Path where initial strain  output is produced


filename = 'wf_junkrad.txt';
delimiter=' ';
data = importdata(filename, delimiter);
for i =218:length(data)	

    wfdata = data(i);
    line = regexp(wfdata, '\t', 'split');
    chararr = line{1};
    gtid = chararr(1);
    wfname = chararr(2);
    
    Y = sc_search('keyword',wfname, 'verbose');
    wfinfo = Y.init_info;
    [~,name,ext] = fileparts(wfinfo.par_file_string);
    wfpath = wfinfo.sim_dir;
    fprintf('%d. Creating Strain for %s waveform - %s \n \n ',i-1, name, gtid{:});
    %try
    func_Strain(outdirpath, wfpath, numrelpath);
   % catch
   %     clearvars -except i, data, outdirpath, wfpath, numrelpath
   %     kit_dir    = '~/DAKIT/'; % location of your analysis toolkit
   %     path(genpath(kit_dir),path);       % add all directories in kit_dir to path

   %     healykit_dir    = '~/Documents/MATLAB/HealyToolkit/'; % location of your analysis toolkit
   %     path(genpath(healykit_dir),path);       % add all directories in kit_dir to path

   % 	i = i-1 ;
   % end
    
 end
%Waveform with Issues -  
%(90) Sq2_d6.2_a0.6_oth.030_rr_M140 (GT0431)
%(100) Sq2_d6.2_a0.6_oth.270_rr_M140
%(217)  D7.5_q15.00_a0.0_CHgEEB_m800 waveform  (GT0601)
%(404) SO_D9_q1.5_th2_135_ph1_0_m120 waveform (GT0847)
%(412) SO_D9_q1.5_th2_135_ph1_135_m140 (GT0858)
%(425) SO_D9_q1.5_th2_135_ph1_180_m140 (GT0871)
%(414)  D8_q1.0_a0.6_0.6_th90_0_r100_res140_CE (GT0860)
%(415) D8_q1.0_a0.6_0.6_th90_135_r100_res140_CE (GT0861)
%(424) SO_D9_q1.5_th2_135_ph1_90_m140 (GT0870)
%(425) SO_D9_q1.5_th2_135_ph1_90_m140 (GT0871) - Problem with interpolation due to duplicate time entries in psi4_analysis file. Either select the unique values for time and remove the corresponding values from other quantities in bn.rundata.system or recompile the simulation and remove the duplicate time values from psi4analysis files. 
%(444) D9_q2.0_a0.6_0.6_th0_r100_res140_CE (GT0890)
%(445) D9_q2.0_a0.6_0.6_th135_r100_res140_CE (GT0891)
%(447) D9_q2.0_a0.6_0.6_th270_r100_res140_CE (GT0894)
%(449) D9_q2.0_a0.6_0.6_th45_r100_res140_CE (GT0896)
%(450) D9_q2.0_a0.6_0.6_th90_r100_res140_CE (GT0897)
