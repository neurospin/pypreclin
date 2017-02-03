% loop trough files in folder to make undistortion for all .nii files in
% folder

% Attention: the mask needs to be made first! (see
% Use_mc_undist_wouter_dante.m)
% This is STEP 2.

clear
% addpaths
addpath('/fmri/apps/spm5')
addpath('/fmri/spm5_scripts/')
addpath('/fmri/spm5_templates/')
addpath('/fmri/spm5_utils/')
% 1- setpath
% for preprocessing
% addpath('/fmri/spm5_utils/prepro_v4')


addpath('/fmri/spm5_utils/prepro_v5_john/')

% 1- setpath
addpath('/fmri/spm5_utils/prepro_v5_john/prepro_tools/prod/mex')
addpath('/fmri/spm5_utils/prepro_v5_john/prepro_tools/prod/undist')

% set the dir where the analysis files are: (must only contain .nii files
% (with 'run' in the name), which have to be u.distorted!)
%  edir= '/data/fmri_monkey_03/PROJECT/John/HighRes_Functional/rawdata/Nancy120803/divideTest/'
%  cd(edir)
%  
%  mc_undist_wouter_dante_natalie('cr_Nancy120803_run000*','cr_Nancy120803_gre3d_resamp_brain_mask.nii');


edir= '/data/fmri_monkey_03/PROJECT/John/HighRes_Functional/rawdata/Stitch121027/'
dd2=dir([edir 'cr*_run*.nii'])
runn=13:16;

%mc_undist_wouter_dante(edir,runn);
mc_undist_wouter_dante_natalie(edir,runn);