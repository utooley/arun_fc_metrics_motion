%% THIS SCRIPT RUNS MODULARITY ANALYSIS ON 4 RUNS of EACH SUBJECT USING GORDON PARCELLATION
% FIX MATRICES ONLY

%% SETUP
%mounted locally on work computer
subj_dir='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/Covariates'
data_dir='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/FunctionalConnectivityMatrices'
scripts_dir=
outdir='~/Documents/projects/in_progress/arun_fc_metrics_motion/output/data/'

%subject list
subjList=readtable(fullfile(subj_dir, 'S1200_Release_Subjects_Demographics.csv'));
subjList=subjList.Subject;

rest_runs={'_REST1_LR_', '_REST1_RL_','_REST2_LR_','_REST2_RL_'};
fc_metrics={'Coherence', 'MutualInformation', 'MutualInformationTime','Pearson','Spearman', 'WaveletCoherence'};
%% initialize vectors
num_communities=zeros(length(subjList),1);

%%%%%%
%for each of four runs
%%%%%%
for i=1:4
    run=rest_runs{i};
%%%%%%
%for each FC metric
%%%%%%
for j=1:6
    metric=fc_metrics{j}
 
% while the average number of communities detected across participants is
% less than 7, keep iterating
avgnumcommunities=0
    while (avgnumcommunities < 6.5 | avgnumcommunities > 7.5)
        gamma=2;
%%%%%%
%for each subject
%%%%%%
for n=1:length(subjList);
    sub=subjList(n)
    try
    file=fullfile(data_dir,strcat('gordon_',num2str(sub),run,'FIX_matrices_', metric,'.mat'));
    load(file);
    avgweight.(metric)(n,1)=mean(AdjMat(AdjMat~=0)); %get average weight for each metric
%% calculate the modularity quality index raw on each metric
%probably won't use this, but worth having
if (j==4 | j == 5) %Pearson or spearman, use the negative weighting  
    %Community Louvain outputs a measure of modularity and can take signed
    %networks as input. Weighted the negative connections asymmetrically, Q* as
    %recommended by Rubinov & Sporns
    [M Q]=community_louvain(AdjMat, gamma, [], 'negative_asym');
    modul(n,1)=Q;
    num_communities(n,1)=length(unique(M)); %how many communities were output
else
    %use the default modularity
    [M Q]=community_louvain(AdjMat, gamma);
    modul(n,1)=Q;
    num_communities(n,1)=length(unique(M));
end
    catch
    disp('This subject not found')
    end
end
avgnumcommunities=mean(num_communities(num_communities~=0))
gamma=gamma+0.3
    end
end

%% Save outfiles for each run
outfile=dataset(avgweight.Pearson, avgweight.Spearman, avgweight.Coherence, avgweight.WaveletCoherence, avgweight.MutualInformation, avgweight.MutualInformationTime, modul, num_communities, gamma)
filename=strcat('modularity_raw',run, '51619.csv')
export(outfile,'File',fullfile(outdir,filename),'Delimiter',',')
    end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Reranking edge weights and calculate modularity %%%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
%for each of four runs
%%%%%%
for i=1:4
    run=rest_runs{i};
    
%%%%%%
%for each FC metric
%%%%%%
for j=1:6
metric=fc_metrics{j}
num_communities.(metric)=zeros(length(subjList),1); %set up num communities

% while the average number of communities detected across participants is
% less than 7, keep iterating
    avgnumcommunities=0
    gamma=1;
    while (avgnumcommunities < 6.5 | avgnumcommunities > 7.5)
%%%%%%
%for each subject
%%%%%%
for n=1:length(subjList);
    sub=subjList(n)
    try
    %load the Pearson correlation matrix
    file=fullfile(data_dir,strcat('gordon_',num2str(sub),run,'FIX_matrices_Pearson.mat'));
    load(file);
    AdjMat( AdjMat == 0 ) = NaN; ; %to avoid 0's being sorted in the order
    PearsMat=AdjMat;
    %pull out the weights in the Pearson correlation (or whatever base)
    %matrix in order from highest to lowest
    num_edges=(size(PearsMat,2))*(size(PearsMat,2)-1)/2 ;%how many?
    weights= maxk(PearsMat(:),(num_edges*2)); %get them, each edge is twice in the adjmat, they are 1=biggest to smallest
    weights_unique= weights(1:2:end);%remove every other row
    %replace the rank ordered weights in the xxx correlation matrix with those
    %Pearson correlation matrix weights
    %load the matrix
    file=fullfile(data_dir,strcat('gordon_',num2str(sub),run,'FIX_matrices_', metric,'.mat'));
    load(file);
    %AdjMat( AdjMat == 0 ) = NaN; %to avoid 0's being sorted in the order
    values_only=sort(AdjMat(AdjMat~=0),'descend'); %no 0's being sorted in the order
    values_only_unique=values_only(1:2:end);%remove duplicate rows from each edge appearing twice
    [l,idx]=ismember(AdjMat,values_only_unique); %order the edge weights by weight, index is 1=largest to smallest
    AdjMat_swap=AdjMat; %make a copy of the matrix
    size(unique(idx(:)))%check length of the index, should match weights_unique!
    for i=1:length(weights_unique)    %put them back in AdjMat_swap in the order of idx, put the Pearson weights into this matrix, in order
        AdjMat_swap(idx==i)=weights_unique(i);
    end
    %make sure the diagonal is zero just in case!
    for x=1:333
        AdjMat_swap(x,x)=0;
    end
    %calculate modularity on the re-weighted xxx matrix, using assymmetric
    %weighting of negative weights since Pearson
    [M Q]=community_louvain(AdjMat_swap, gamma, [], 'negative_asym');
    modul.(metric)(n,1)=Q;
    num_communities.(metric)(n,1)=length(unique(M));
    num_communities.(metric)%how many communities were output
    %cycle through to ensure that the mean number of communities is 7
    %(mean across iterations of Gen Louvain, or mean across subjects)?
    
    %%save outfile for each metric
    %WHY DOES PEARSON AND SPEARMAN HAVE MORE UNIQUE WEIGHTS THAN EDGES?
%     %Troubleshooting
%     [C, ia, ic] = unique(PearsMat,'stable')
%     dups=reshape(ic, size(PearsMat))

  %% CHECKS FOR MATRIX SWAPS
%     find(AdjMat_swap == max(AdjMat_swap(:)))
%     find(AdjMat == max(AdjMat(:)))
%     find(AdjMat_swap == min(AdjMat_swap(AdjMat_swap~=0)))
%     find(AdjMat == min(AdjMat(AdjMat~=0))) %these work, great! Except for Pearson matrices, which have overlapping unique values > num_edges. 
    %Figured out that it was because of floating point precision of the
    %Pearson and Spearman matrices
    catch
        disp('This subject not found')
    end
    end
    avgnumcommunities=mean(num_communities.(metric)(num_communities.(metric)~=0))
    gamma=gamma+0.3
end
end
    
%% Save outfiles for each run
save(fullfile(outdir, 'modul'), modul.(metric))
save(fullfile(outdir, 'numcommunities'), num_communities.(metric))

outfile=dataset(modul.Pearson, modul.Spearman, modul.Coherence, modul.WaveletCoherence, modul.MutualInformation, modul.MutualInformationTime, num_communities.Pearson, num_communities.Spearman, num_communities.Coherence, num_communities.WaveletCoherence, num_communities.MutualInformation, num_communities.MutualInformationTime, gamma)
filename=strcat('modularity_reranked_Pearson_weights_into',metric, run, '51619.csv')
export(outfile,'File',fullfile(outdir,filename),'Delimiter',',')


end

