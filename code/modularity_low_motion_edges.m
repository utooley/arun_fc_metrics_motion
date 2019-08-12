%% THIS SCRIPT RUNS MODULARITY ANALYSIS ON 4 RUNS of EACH SUBJECT USING GORDON PARCELLATION
% FIX MATRICES ONLY

%% SETUP
%mounted locally on work computer
subj_dir='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/Covariates'
data_dir='~/Desktop/cluster/jag/bassett-lab/hcp_Max/Data/QCFC_correlationMatrices'
outdir='~/Documents/projects/in_progress/arun_fc_metrics_motion/output/data/'
%working directly on the cluster
subj_dir='/data/jag/bassett-lab/hcp_Max/Data/Covariates'
data_dir='/data/jag/bassett-lab/hcp_Max/Data/QCFC_correlationMatrices'
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/arun_fc_metrics_motion/output/data/Gordon_ICA_FIX'
%working locally on personal computer
subj_dir='~/Documents/projects/in_progress/arun_fc_metrics_motion/data/subjLists/'
data_dir='~/Documents/projects/in_progress/arun_fc_metrics_motion/data/QCFC_correlationMatrices'
outdir='~/Documents/projects/in_progress/arun_fc_metrics_motion/output/data/'
%subject list
subjList=readtable(fullfile(subj_dir, 'S1200_Release_Subjects_Demographics.csv'));
subjList=subjList.Subject;

rest_runs={'_REST1_LR_', '_REST1_RL_','_REST2_LR_','_REST2_RL_'};
fc_metrics={'Coherence', 'MutualInformation', 'MutualInformationTime','Pearson','Spearman', 'WaveletCoherence'};
%% initialize vectors

%%%%%
%for each of four runs
%%%%%%
%find the edges across all 6 metrics that are the least overall affected by
%motion
for i=1:4
    matrix=zeros(333,333);
    run=rest_runs{i};
    run_name=strcat('run',run, 'l');
for j=1:6 %only Pearson and Spearman
    metric=fc_metrics{j}
    %read in the average matrix, take the absolute value of correlation,
    %and rank it
    file=fullfile(data_dir,strcat('QCFC_correlationMatrix_',metric,'_gordon',run,'FIX_matrices.mat'));
    load(file);
    QCFC_correlationMatrix_abs=abs(QCFC_correlationMatrix);
    [~, ranks]=ismember(QCFC_correlationMatrix_abs,unique(sort(QCFC_correlationMatrix_abs,'descend')));
    %% FIX THIS, THE RANKING ISN"T WORKING. MAYBE NOW?
    matrix=matrix+ranks; %add each ranked matrix to the other.
end
%threshold matrix to the 20% of edges with the smallest values
thresholded=unique(sort(matrix(:),'ascend'))       
[~, ranks]=ismember(matrix,thresholded(1:ceil(length(thresholded)*0.2)));
%I think this is working
   
%save out this matrix as a mask for each run
save(fullfile(outdir, strcat('noMotion_edges_','gordon', run,'FIX_matrices.mat'), 'ranks')
end
%%%%%%
%for each of four runs
%%%%%%
for i=1:4
    run=rest_runs{i};
    run_name=strcat('run',run, 'l');
    modul.(metric)=zeros(length(subjList),1);
    avgweight.(metric)=zeros(length(subjList),1);
%%%%%%
%for each FC metric
%%%%%%
for j=4:5 %only Pearson and Spearman
    metric=fc_metrics{j}
 num_communities.(metric)=zeros(length(subjList),1); %set up num communities
 
% while the average number of communities detected across participants is
% less than 7, keep iterating
avgnumcommunities=0
gamma=1;
    while (avgnumcommunities < 6 | avgnumcommunities > 8)
%%%%%%
%for each subject
%%%%%%
for n=1:length(subjList);
    sub=subjList(n);
    try
    file=fullfile(data_dir,strcat('gordon_',num2str(sub),run,'FIX_matrices_', metric,'.mat'));
    load(file);
    AdjMat=threshold_absolute(AdjMat,0);
    avgweight.(metric)(n,1)=mean(AdjMat(AdjMat~=0)); %get average weight for each metric
%% calculate the modularity quality index raw on each metric
%probably won't use this, but worth having
% if (j==4 | j == 5) %Pearson or spearman, use the negative weighting  
%     %Community Louvain outputs a measure of modularity and can take signed
%     %networks as input. Weighted the negative connections asymmetrically, Q* as
%     %recommended by Rubinov & Sporns
%     [M Q]=community_louvain(AdjMat, gamma, [], 'negative_asym');
%     modul.(metric)(n,1)=Q;
%     num_communities.(metric)(n,1)=length(unique(M)); %how many communities were output
% else
    %use the default modularity
    [M Q]=community_louvain(AdjMat, gamma);
    modul.(metric)(n,1)=Q;
    num_communities.(metric)(n,1)=length(unique(M));
% end
    catch
    disp('This subject not found')
    end
end
avgnumcommunities=mean(num_communities.(metric)(num_communities.(metric)~=0))
gamma=gamma+0.01
    end
    allgamma.(metric).(run_name)=gamma
end
%% Save outfiles for each run
% outfile=dataset(avgweight.Pearson, avgweight.Spearman, avgweight.Coherence, avgweight.WaveletCoherence, avgweight.MutualInformation, avgweight.MutualInformationTime, modul.Pearson, modul.Spearman, modul.Coherence, modul.WaveletCoherence, modul.MutualInformation, modul.MutualInformationTime, num_communities.Pearson, num_communities.Spearman, num_communities.Coherence, num_communities.WaveletCoherence, num_communities.MutualInformation, num_communities.MutualInformationTime,  repmat(allgamma.Pearson.(run_name),length(subjList),1), repmat(allgamma.Spearman.(run_name),length(subjList),1), repmat(allgamma.Coherence.(run_name),length(subjList),1), repmat(allgamma.WaveletCoherence.(run_name),length(subjList), 1), repmat(allgamma.MutualInformation.(run_name), length(subjList),1), repmat(allgamma.MutualInformationTime.(run_name), length(subjList), 1))
% header={'avgweight_Pearson',	'avgweight_Spearman',	'avgweight_Coherence',	'avgweight_WaveletCoherence',	'avgweight_MutualInformation',	'avgweight_MutualInformationTime',	'modul_Pearson',	'modul_Spearman',	'modul_Coherence',	'modul_WaveletCoherence',	'modul_MutualInformation',	'modul_MutualInformationTime'	,'num_communities_Pearson'	,'num_communities_Spearman'	,'num_communities_Coherence'	,'num_communities_WaveletCoherence','num_communities_MutualInformation'	,'num_communities_MutualInformationTime'	,'allgamma_Pearson'	,'allgamma_Spearman',	'allgamma_Coherence','allgamma_WaveletCoherence','allgamma_MutualInformation','allgamma_MutualInformationTime'}
% outfile.Properties.VarNames=header
% filename=strcat('modularity_raw',run, '071619.csv') %need to figure out how to get the headers to work.
% export(outfile,'File',fullfile(outdir,filename),'Delimiter',',')

 
save(fullfile(outdir, strcat('modul_positiveonly_run',int2str(i))), 'modul')
save(fullfile(outdir, strcat('numcommunities_positiveonly_run',int2str(i))), 'num_communities')
save(fullfile(outdir, strcat('avgnumcommunities_positiveonly_run',int2str(i))), 'avgnumcommunities')
save(fullfile(outdir, strcat('gamma_positiveonly_run',int2str(i))), 'allgamma')
    
end
%save the all gamma variables, just in case
filename=fullfile(outdir,'allgamma_modularity_raw_071619.mat')
save(filename, 'allgamma')



%%%%%%%%%%%%%%%%%%%%%%%
%%% Reranking edge weights and calculate modularity %%%
%%%%%%%%%%%%%%%%%%%%%%
clear num_communities
baseline_metric={'Pearson','WaveletCoherence'};
fc_metrics={'Pearson','WaveletCoherence'};

%may need to reload allgamma .mat file here
load(fullfile(outdir,'allgamma_modularity_raw.mat'))
%%%%%%
%for each of four runs
%%%%%%
for i=1:4
    run=rest_runs{i};
    run_name=strcat('run',run, 'l');
    clear modul
%%%%%%
%for each subject
%%%%%%
for n=1:length(subjList);
    sub=subjList(n);
    %%%%%%
    %for each of two baseline edge distributions
    %%%%%%
    for l=1:2
        baseline_met=baseline_metric{l};
        % while the average number of communities detected across participants is
        % less than 7, keep iterating
        %avgnumcommunities=0
        %gamma=1;
        %while (avgnumcommunities < 6.5 | avgnumcommunities > 7.5)
        %or just pull the gamma from above, from the tuning for the edge
        %distribution being the one above
        gamma=allgamma.(baseline_met).(run_name); %SET GAMMA FROM BEFORE
    %%%%%%
    %for each FC metric
    %%%%%%
    for j=1:2
    metric=fc_metrics{j}
    try
    %load the Pearson correlation matrix
    file=fullfile(data_dir,strcat('gordon_',num2str(sub),run,'FIX_matrices_',baseline_met,'.mat'));
    load(file);
    AdjMat( AdjMat == 0 ) = NaN; %to avoid 0's being sorted in the order
    PearsMat=AdjMat;
    %pull out the weights in the Pearson correlation (or whatever base)
    %matrix in order from highest to lowest
    num_edges=(size(PearsMat,2))*(size(PearsMat,2)-1)/2 ;%how many?
    %minmax_install needed, ensure code folder is on matlab path
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
    [lia,idx]=ismemberf(AdjMat,values_only_unique); %order the edge weights by weight, index is 1=largest to smallest
    %this is a problem for Pearson matrices--some of the values don't match
    %those in values_only_unique because of floating point tolerance. Use
    %ismemberf instead.
    AdjMat_swap=AdjMat; %make a copy of the matrix
    size(unique(idx(:)))%check length of the index, should match weights_unique! Off by 1 because idx has 0's on the diagonal, and weights unique has no 0s.
    for w=1:length(weights_unique)    %put them back in AdjMat_swap in the order of idx, put the Pearson weights into this matrix, in order
        AdjMat_swap(idx==w)=weights_unique(w);
    end
    %make sure the diagonal is zero just in case!
    for x=1:333
        AdjMat_swap(x,x)=0;
    end
    %calculate modularity on the re-weighted xxx matrix, using assymmetric
    %weighting of negative weights since Pearson
    if (strcmp(baseline_met,'Pearson'))
        [M Q]=community_louvain(AdjMat_swap, gamma, [], 'negative_asym');
    else
        [M Q]=community_louvain(AdjMat_swap, gamma);
    end
    modul.(baseline_met).(metric)(n,1)=Q;
    num_communities.(baseline_met).(metric)(n,1)=length(unique(M));
    %num_communities.(baseline_met).(metric)%how many communities were output
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
    
    save(fullfile(outdir, strcat('modul_run',int2str(i))), 'modul')
    save(fullfile(outdir, strcat('numcommunities_run',int2str(i))), 'num_communities')
    save(fullfile(outdir, strcat('avgnumcommunities_run',int2str(i))), 'avgnumcommunities')
    catch
        disp('This subject not found')
    end
    end
    
    end
end
    avgnumcommunities.(baseline_met).(metric)=mean(num_communities.(baseline_met).(metric)(num_communities.(baseline_met).(metric)~=0))

    
%% Save outfiles
try
outfile=dataset(modul.Pearson.Pearson, modul.Pearson.Spearman,modul.Pearson.Coherence, modul.Pearson.WaveletCoherence, modul.Pearson.MutualInformation, modul.Pearson.MutualInformationTime, modul.WaveletCoherence.Pearson, modul.WaveletCoherence.Spearman,modul.WaveletCoherence.Coherence, modul.WaveletCoherence.WaveletCoherence, modul.WaveletCoherence.MutualInformation, modul.WaveletCoherence.MutualInformationTime,num_communities.Pearson.Pearson, num_communities.Pearson.Spearman, num_communities.Pearson.Coherence, num_communities.Pearson.WaveletCoherence, num_communities.Pearson.MutualInformation, num_communities.Pearson.MutualInformationTime,num_communities.WaveletCoherence.Pearson, num_communities.WaveletCoherence.Spearman, num_communities.WaveletCoherence.Coherence, num_communities.WaveletCoherence.WaveletCoherence, num_communities.WaveletCoherence.MutualInformation,num_communities.WaveletCoherence.MutualInformationTime)

%fix the header here
header={'avgweight_Pearson',	'avgweight_Spearman',	'avgweight_Coherence',	'avgweight_WaveletCoherence',	'avgweight_MutualInformation',	'avgweight_MutualInformationTime',	'modul_Pearson',	'modul_Spearman',	'modul_Coherence',	'modul_WaveletCoherence',	'modul_MutualInformation',	'modul_MutualInformationTime'	,'num_communities_Pearson'	,'num_communities_Spearman'	,'num_communities_Coherence'	,'num_communities_WaveletCoherence','num_communities_MutualInformation'	,'num_communities_MutualInformationTime'	,'allgamma_Pearson'	,'allgamma_Spearman',	'allgamma_Coherence','allgamma_WaveletCoherence','allgamma_MutualInformation','allgamma_MutualInformationTime'}
outfile.Properties.VarNames=header

filename=strcat('modularity_reranked_Pearson_weights_into',metric, run, '7219.csv')
export(outfile,'File',fullfile(outdir,filename),'Delimiter',',')

%export average number of communities
outfile2=table(avgnumcommunities.Pearson.Pearson, avgnumcommunities.Pearson.Spearman,av.numcommunities.Pearson.Coherence, avgnumcommunities.Pearson.WaveletCoherence, avgnumcommunities.Pearson.MutualInformation, avgnumcommunities.Pearson.MutualInformationTime, avgnumcommunities.WaveletCoherence.Pearson, avgnumcommunities.WaveletCoherence.Spearman, avgnumcommunities.WaveletCoherence.Coherence, avgnumcommunities.WaveletCoherence.WaveletCoherence, avgnumcommunities.WaveletCoherence.MutualInformation, avgnumcommunities.WaveletCoherence.MutualInformationTime)
save(fullfile(outdir, 'avgnumcommunities_071619'),'outfile2')
catch
    
    disp('saving csvs didnt work')
end
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate modularity for each metric with respect to the a prior community partition %%%
%%%%%%%%%%%%%%%%%%%%%%
clear modul
clear modul2
clear numcommunities
%For Gordon
%sparcels=readtable('~/Documents/toolboxes/parcellations/Gordon_333/Parcels.csv')
parcels=readtable('/data/picsl/mackey_group/tools/gordon333/Parcels.csv')
parcels=parcels(:,{'Community','ParcelID'})
parcels.Community=categorical(parcels.Community)
%For Yeo
parcels=readtable('~/Documents/toolboxes/parcellations/')

%%%%%%
%for each of four runs
%%%%%%
for i=1:4
    run=rest_runs{i};
    run_name=strcat('run',run, 'l');
%%%%%%
%for each FC metric
%%%%%%
for j=1:6
    metric=fc_metrics{j}
    modul.(metric)=zeros(length(subjList),1);
    modul2.(metric)=zeros(length(subjList),1);
%%%%%%
%for each subject
%%%%%%
    for n=1:length(subjList);
        sub=subjList(n);
        try
        file=fullfile(data_dir,strcat('gordon_',num2str(sub),run,'FIX_matrices_', metric,'.mat'));
        load(file);
        AdjMat=threshold_absolute(AdjMat,0);
        avgweight.(metric)(n,1)=mean(AdjMat(AdjMat~=0)); %get average weight for each metric
    %% calculate the modularity quality index raw on each metric
    %probably won't use this, but worth having
%     if (j==4 | j == 5) %Pearson or spearman, use the negative weighting  
%         %Community Louvain outputs a measure of modularity and can take signed
%         %networks as input. Weighted the negative connections asymmetrically, Q* as
%         %recommended by Rubinov & Sporns
%         %[M Q]=community_louvain(AdjMat, [], double(parcels.Community), 'negative_asym');
%         [M Q]=modul_only(AdjMat, [], double(parcels.Community), 'negative_asym');
%         modul.(metric)(n,1)=Q;
%         Q %check that it's still the same
%         num_communities.(metric)(n,1)=length(unique(M)); %how many communities were output
%     else
        %use the default modularity
        %[M Q]=community_louvain(AdjMat, [], double(parcels.Community));
        [M Q]=modul_only(AdjMat, [], double(parcels.Community));
        modul.(metric)(n,1)=Q; %check that it's still the same.
        num_communities.(metric)(n,1)=length(unique(M));
        Q =QFModul(double(parcels.Community), AdjMat); %try a different modularity implementation
        modul2.(metric)(n,1)=Q;
%     end
        catch
        disp('This subject not found')
        end
    end
end
%% Save outfiles for each run
save(fullfile(outdir, strcat('modul_pose_gordon_partition',int2str(i))), 'modul')
save(fullfile(outdir, strcat('modul2_pos_gordon_partition',int2str(i))), 'modul2')
save(fullfile(outdir, strcat('numcomms_pos_gordon_partition',int2str(i))), 'num_communities')
try
    outfile=dataset(avgweight.Pearson, avgweight.Spearman, avgweight.Coherence, avgweight.WaveletCoherence, avgweight.MutualInformation, avgweight.MutualInformationTime, modul.Pearson, modul.Spearman, modul.Coherence, modul.WaveletCoherence, modul.MutualInformation, modul.MutualInformationTime, num_communities.Pearson, num_communities.Spearman, num_communities.Coherence, num_communities.WaveletCoherence, num_communities.MutualInformation, num_communities.MutualInformationTime,  repmat(allgamma.Pearson.(run_name),length(subjList),1), repmat(allgamma.Spearman.(run_name),length(subjList),1), repmat(allgamma.Coherence.(run_name),length(subjList),1), repmat(allgamma.WaveletCoherence.(run_name),length(subjList), 1), repmat(allgamma.MutualInformation.(run_name), length(subjList),1), repmat(allgamma.MutualInformationTime.(run_name), length(subjList), 1))
    header={'avgweight_Pearson',	'avgweight_Spearman',	'avgweight_Coherence',	'avgweight_WaveletCoherence',	'avgweight_MutualInformation',	'avgweight_MutualInformationTime',	'modul_Pearson',	'modul_Spearman',	'modul_Coherence',	'modul_WaveletCoherence',	'modul_MutualInformation',	'modul_MutualInformationTime'	,'num_communities_Pearson'	,'num_communities_Spearman'	,'num_communities_Coherence'	,'num_communities_WaveletCoherence','num_communities_MutualInformation'	,'num_communities_MutualInformationTime'	,'allgamma_Pearson'	,'allgamma_Spearman',	'allgamma_Coherence','allgamma_WaveletCoherence','allgamma_MutualInformation','allgamma_MutualInformationTime'}
    outfile.Properties.VarNames=header
    filename=strcat('modularity_Gordon_partition',run, '072719.csv') %need to figure out how to get the headers to work.
    export(outfile,'File',fullfile(outdir,filename),'Delimiter',',')
catch
    
    disp('saving csvs didnt work')
end
end
