%% SETUP
%mounted locally on work computer
subj_dir='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/Covariates'
data_dir='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/FunctionalConnectivityMatrices'
scripts_dir=

%subject list
subjList=readtable(fullfile(subj_dir, 'S1200_Release_Subjects_Demographics.csv'));
subjList=subjList.Subject;
%%
%for each subject

%for each FC metric

%calculate the modularity quality index raw on each metric
%probably won't use this, but worth having

    %cycle through to ensure that the mean number of communities is 7
    %(mean across iterations of Gen Louvain, or mean across subjects)?

%rank order the weights in the Pearson correlation matrix

%replace the rank ordered weights in the xxx correlation matrix with the
%Pearson correlation matrix

%calculate modularity on the re-weighted xxx matrix