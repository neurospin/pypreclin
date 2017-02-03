addpath {SPMPATH}
addpath {PREPROCPATH}/resources/leuven/prod/undist
addpath {PREPROCPATH}/resources/leuven/prod/mex/

currentDir = '{CWD}';
mc_undist_wouter_dante_tasserie(currentDir, '{BASENAME}')
exit;

