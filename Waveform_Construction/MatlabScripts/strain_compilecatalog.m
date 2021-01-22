%% Script to generate strain data from entire simulation catalog

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
for i =1:length(data)	

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
    func_Strain(outdirpath, wfpath, numrelpath);
    
 end
