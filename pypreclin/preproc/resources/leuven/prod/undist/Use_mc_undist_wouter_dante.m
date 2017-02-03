% use mc_untidst_wouter_dante.m
% ATTENTION: this program is needed to remomve NAN values from the to GRE
% registered MASK:

% STEP1a: Do Coreg.Estimate Reslice in SPM (with
% GRE(Reference),Anat(Source), and Mask(Other)
% STEP1b: before using mc_undist_wouter_dante.m: here.

clear
% addpaths
addpath('/fmri/apps/spm5')
addpath('/fmri/spm5_scripts/')
addpath('/fmri/spm5_templates/')
addpath('/fmri/spm5_utils/')
% 1- setpath
addpath('/fmri/spm5_utils/prepro_v5_natalie/')

% 1- setpath
addpath('/fmri/spm5_utils/prepro_v5_natalie/prepro_tools/prod/undist')
addpath('/fmri/spm5_utils/prepro_v5_natalie/prepro_tools/prod/mex')


% 2- make mask: 
% - take nancy_cropped_centered_1mm.nii, MRIread, and set all
% - above zero to 1 (with find command) -> this will be my mask to use:
% MRIwrite(a,'..msk.nii')

% in SPM, ESTIMATE&RESLICE: coregister anatomy (nancy_cropped_centered_1mm.nii) with gre (.nii, as converted from
% scanner)
% - Reference: gre
% - Source: anatomy cropped centered
% - Other: ...msk.nii (parameters are applied to mask too!)

% attention: this steps puts NaN values into mask volume, set them back to
% zero doing this (original mask with NaN values will be overwritten)

% set my analysis dir:
edir='/data/fmri_monkey_raw/DICOM/Nancy/Nancy120517/';
dd=dir([edir '*_msk*.nii']);

path_mask=[edir dd(1).name];

V_mask=spm_vol(path_mask);

mask=spm_read_vols(V_mask);
% this figure shows the NaN values at the edges of the volume(if present)
figure ; imagesc(mask(:,:,30));
% check how many NaN values in mask!
sum(isnan(mask(:)))
% V_mask = 
%       fname: [1x79 char]
%         mat: [4x4 double]
%         dim: [96 96 48]
%          dt: [16 0]
%       pinfo: [3x1 double]
%           n: [1 1]
%     descrip: 'spm - realigned'
%     private: [1x1 nifti]
V_mask.fname % gives the name and path of mask
%% this command did not work:
% mask(isnan(mask))=0;
% figure ; imagesc(mask(:,:,30));
%% rather: mask= is all above 0.5, rest is set to zero:
mask=mask>0.5;
figure ; imagesc(mask(:,:,30));
spm_write_vol(V_mask,mask);

