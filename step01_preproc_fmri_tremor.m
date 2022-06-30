clear
clc


%% Prepare paths and regexp

maindir = '/network/lustre/iss02/cenir/analyse/irm/users/benoit.beranger/VIBRALCA';

par.redo= 0;
par.run = 1;
par.pct = 0;
par.sge = 0;


%% Get files paths

e = exam(maindir,'nifti','VIBRALCA');


%% Get files paths #matvol

%Anat
e.addSerie('3DT1$'   , 'anat_T1', 1 );
e.addSerie('3DFLAIR$', 'anat_FLAIR', 1 );
e.getSerie('anat').addVolume('^v_.*nii','v',1);

% Func
e.addSerie('Run01$', 'run_01', 1)
e.addSerie('Run02$', 'run_02', 1)

e.getSerie('run').addVolume('^v_.*nii$', 'v', 1);

e.reorderSeries('name');

e.explore


%% Segment anat with cat12

par.subfolder = 0;         % 0 means "do not write in subfolder"
par.biasstr   = 0.5;
par.accstr    = 0.5;
par.GM        = [1 0 1 0]; %                          (wp1*)     /                        (mwp1*)     /              (p1*)     /                            (rp1*)
par.WM        = [1 0 1 0]; %                          (wp2*)     /                        (mwp2*)     /              (p2*)     /                            (rp2*)
par.CSF       = [1 0 1 0]; %                          (wp3*)     /                        (mwp3*)     /              (p3*)     /                            (rp3*)
par.TPMC      = [1 0 1 0]; %                          (wp[456]*) /                        (mwp[456]*) /              (p[456]*) /                            (rp[456]*)
par.label     = [1 0 0] ;  % native (p0*)  / normalize (wp0*)  / dartel (rp0*)       This will create a label map : p0 = (1 x p1) + (3 x p2) + (1 x p3)
par.bias      = [1 1 0] ;  % native (ms*)  / normalize (wms*)  / dartel (rms*)       This will save the bias field corrected  + SANLM (global) T1
par.las       = [0 0 0] ;  % native (mis*) / normalize (wmis*) / dartel (rmis*)      This will save the bias field corrected  + SANLM (local) T1
par.warp      = [1 1];     % Warp fields  : native->template (y_*) / native<-template (iy_*)
par.doSurface = 0;
par.doROI     = 0;         % Will compute the volume in each atlas region
par.jacobian  = 0;         % Write jacobian determinant in normalize space

anat = e.gser('anat_T1').gvol('^v');
job_do_segmentCAT12(anat,par);

par.jobname = 'zipWMCSF';
e.gser('anat_T1').gvol('^wp[23]').zip_and_keep(par);
par = rmfield(par,'jobname');


%% Preprocess fMRI runs

%realign
par.type = 'estimate';
ffunc_nm = e.getSerie('run').getVolume('^v');
j_realign_reslice_nm = job_realign(ffunc_nm,par);

%coregister mean fonc on brain_anat
fanat = e.getSerie('anat_T1').getVolume('^p0');
fmean = e.getSerie('run').getVolume('^meanv'); fmean = fmean(:,1); % use the mean of the run1 to estimate the coreg
fo    = e.getSerie('run').getVolume('^v');
par.type = 'estimate';
par.jobname = 'spm_coreg_epi2anat';
j_coregister=job_coregister(fmean,fanat,fo,par);
par = rmfield(par,'jobname');

%apply normalize
par.vox    = [2 2 2];
to_warp    = fo.removeEmpty + fmean.removeEmpty;
warp_field = to_warp.getExam.getSerie('anat_T1').getVolume('^y');
j_apply_normalize=job_apply_normalize(warp_field,to_warp,par);

%smooth the data
ffonc = e.getSerie('run').getVolume('wv').removeEmpty;
par.smooth = [4 4 4];
j_smooth=job_smooth(ffonc,par);

% coregister WM & CSF on functionnal (using the warped mean)
if isfield(par,'prefix'), par = rmfield(par,'prefix'); end
ref = e.getSerie('run');
ref = ref(:,1).getVolume('^wmeanv'); % first acquired run (time)
src = e.getSerie('anat_T1').getVolume('^wp2');
oth = e.getSerie('anat_T1').getVolume('^wp3');
par.type = 'estimate_and_write';
par.jobname = 'spm_coreg_wmcsf2epi';
job_coregister(src,ref,oth,par);
par = rmfield(par,'jobname');

% robust EPI mask for firstlevel
epi = e.getSerie('run').getVolume('^wv').removeEmpty;
par.fsl_output_format = 'NIFTI';
do_fsl_robust_mask_epi(epi,par);

save('e','e')


%% PhysIO nuisance regressor generation #matlab/TAPAS-PhysIO
%% Prepare files

run  = e.getSerie('run').removeEmpty;

volume = run.getVolume('^wv');

outdir = volume.getDir();

rp = run.getRP('rp_spm');

mask = run.getExam.getSerie('anat').getVolume('^rwp[23]');


%% Prepare job

%----------------------------------------------------------------------------------------------------------------------------------------------------
% ALWAYS MANDATORY
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.physio   = 0;
par.noiseROI = 1;
par.rp       = 1;

par.TR     = 1.525;
par.nSlice = 58;

par.volume = volume;
par.outdir = outdir;

%----------------------------------------------------------------------------------------------------------------------------------------------------
% Physio
%----------------------------------------------------------------------------------------------------------------------------------------------------

% par.physio_Info = info;
% par.physio_PULS = puls;
% par.physio_RESP = resp;

par.physio_RETROICOR        = 0;
par.physio_HRV              = 0;
par.physio_RVT              = 0;
par.physio_logfiles_vendor  = 'Siemens_Tics'; % Siemens CMRR multiband sequence, only this one is coded yet
par.physio_logfiles_align_scan = 'last';         % 'last' / 'first'
% Determines which scan shall be aligned to which part of the logfile.
% Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans;
par.physio_slice_to_realign    = 'middle';       % 'first' / 'middle' / 'last' / sliceNumber (integer)
% Slice to which regressors are temporally aligned. Typically the slice where your most important activation is expected.


%----------------------------------------------------------------------------------------------------------------------------------------------------
% noiseROI
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.noiseROI_mask   = mask;
par.noiseROI_volume = volume;

par.noiseROI_thresholds   = [0.95 0.80];     % keep voxels with tissu probabilty >= 95%
par.noiseROI_n_voxel_crop = [   2    1];     % crop n voxels in each direction, to avoid partial volume
par.noiseROI_n_components = 10;              % keep n PCA componenets


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Realignment Parameters
%----------------------------------------------------------------------------------------------------------------------------------------------------

par.rp_file = rp;

par.rp_order     = 24;   % can be 6, 12, 24
% 6 = just add rp, 12 = also adds first order derivatives, 24 = also adds first + second order derivatives
par.rp_method    = 'FD'; % 'MAXVAL' / 'FD' / 'DVARS'
par.rp_threshold = 0.5;  % Threshold above which a stick regressor is created for corresponding volume of exceeding value


%----------------------------------------------------------------------------------------------------------------------------------------------------
% Other
%----------------------------------------------------------------------------------------------------------------------------------------------------
par.print_figures = 0; % 0 , 1 , 2 , 3

% cluster
par.jobname  = 'spm_physio';
par.walltime = '04:00:00';
par.mem      = '4G';

job_physio_tapas( par );


save e e
