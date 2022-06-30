clear
clc

load e 


%% Specify

model_name = 'model_1';

dirStats = e.mkdir('glm',model_name);
run = e.getSerie('run');
rp = run.getRP('multiple_regressors');
dirFunc  = get_parent_path( run.getVolume('swv').toJob() );
dirFunc  = cellfun(@cellstr, dirFunc, 'UniformOutput', 0);

onsetspath = '/network/lustre/iss02/cenir/analyse/irm/users/benoit.beranger/VIBRALCA/onsets';


e.getSerie('run_01').addStim(onsetspath, 'run01_SPM.mat$', 'run_01')
e.getSerie('run_02').addStim(onsetspath, 'run02_SPM.mat$', 'run_02')

onsets = e.getSerie('run').getStim().toJob(0);



clear par
par.file_reg = '^swv_.*nii';
par.rp       = 1;
par.rp_regex = '^multiple_regressors.txt';

% Masking
par.mask_thr = 0.1; % spm default option
run1 = run(:,1);
mask_run1 = run1.getVolume('mask').toJob(0);
par.mask     =  mask_run1; % cell(char) of the path for the mask of EACH model : N models means N paths

par.sge      = 0;
par.run      = 1;
par.display  = 0;
par.redo     = 0;

par.TR = 1.580;

% Run job
job_first_level_specify(dirFunc,dirStats,onsets,par);
e.addModel('glm',model_name,model_name);
model_obj = e.getModel(model_name);
fspm = model_obj.getPath();


%% Estimate

clear par
par.write_residuals = 0;

par.jobname  = 'spm_glm_est';
par.walltime = '11:00:00';

par.sge      = 0;
par.run      = 1;
par.display  = 0;
par.redo     = 0;

job_first_level_estimate(fspm, par);


%% Contrast

clear par
par.sessrep         = 'both';
par.report          = 0;

par.jobname         ='spm_glm_con';
par.walltime        = '04:00:00';

par.sge             = 0;
par.run             = 1;
par.display         = 0;
par.delete_previous = 1;


rest  = [1 0 0 0 0];
extD  = [0 1 0 0 0];
flexD = [0 0 1 0 0];
extG  = [0 0 0 1 0];
flexG = [0 0 0 0 1];


contrast_T.values = {
    
extD+flexD - 2*rest
extG+flexG - 2*rest
extD+flexD+extG+flexG - 4*rest

+((extD+flexD) - (extG+flexG))
-((extD+flexD) - (extG+flexG))

extD-rest
flexD-rest
extG-rest
flexG-rest

+(extD-flexD)
-(extD-flexD)
+(extG-flexG)
-(extG-flexG)

3*extD  - (flexD + flexG + extG)
3*flexD - ( extD + flexG + extG)
3*extG  - (flexG + flexD + extD)
3*flexG - ( extG + flexD + extD)

}';

contrast_T.names = {
    
'extD+flexD - 2*rest'
'extG+flexG - 2*rest'
'extD+flexD+extG+flexG - 4*rest'

'+((extD+flexD) - (extG+flexG))'
'-((extD+flexD) - (extG+flexG))'

'extD-rest'
'flexD-rest'
'extG-rest'
'flexG-rest'

'+(extD-flexD)'
'-(extD-flexD)'
'+(extG-flexG)'
'-(extG-flexG)'

'3*extD  - (flexD + flexG + extG)'
'3*flexD - ( extD + flexG + extG)'
'3*extG  - (flexG + flexD + extD)'
'3*flexG - ( extG + flexD + extD)'

}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));

contrast_F.names  = {'F-all'};
contrast_F.values = {eye(5)};
contrast_F.types  = cat(1,repmat({'F'},[1 length(contrast_F.names)]));

contrast.names  = [contrast_F.names  contrast_T.names ];
contrast.values = [contrast_F.values contrast_T.values];
contrast.types  = [contrast_F.types  contrast_T.types ];

job_first_level_contrast(fspm,contrast,par);


%% Show

model_obj(1).show()
