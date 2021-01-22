%% Script to generate strain data for single simulation

clc; clear;

% add HealyToolkit path to matlab search directory

healykit_dir    = '~/Documents/MATLAB/HealyToolkit/'; % location of your analysis toolkit
path(genpath(healykit_dir),path);       % add all directories in healykit_dir to path


%Catalog Path where final waveforms files will be created and stored.
numrelpath = '/numrel/NumRel/bkhamesra3/Finalized_Waveforms/Waveform_files/Remaining/';  %Path on numrel where strain output is moved

%Create a temporary directory and add the path here - the strain files will
%be generated here and then will be moved to numrelpath
outdirpath = '/localdata2/bkhamesra3/Waveforms'; %Path where initial strain  output is produced

% Add the simulation directory which has combined data from all the outputs. 
wfpath = '/nethome/numrel/datafiles/Waveforms/SO-series/SO_D9_q1.5_th2_135_ph1_180_m140';
func_Strain(outdirpath, wfpath, numrelpath);
