% Federated Tensor Factorization for Computational Phenotyping. ACMKDD 17
% Yejin Kim (yejin.kim89@gmail.com)
% 


%% load data
clear;
addpath(genpath('./tensor_toolbox'));
rank=10; %rank
maxiter=100;
repeat=1;
max_count=3;
dataset='count_mimic3';

% Read MIMICIII dataset
fileName=strcat(dataset, '.csv');
count_mat=csvread(fileName); % count, patient, medi, diag
count_mat=uint16(count_mat);

count_mat((count_mat(:,1)>max_count),1)=max_count;

sz=max(count_mat);

eleM=1:sz(3);
eleD=1:sz(4);
numP=sz(2);
guide=zeros(sz(3),1);

%% Configurate num of partitions or skewness. (Default is 2 partitions with 0 skewness)
%% num of partition
skewness=0;
numK=[2,3,4,5];
for K=numK;
    run_partition;
    %run_partition(5, skewness, count_mat, guide, sz, eleM, eleD, numP, repeat, maxiter, rank) ;
    
end


%% skewness
%K=3;
%skewnessVal=[0.5, 0.7, 0.9];
%for s=1:3
%    skewness=skewnessVal(s);
%    run_partition;
    %run_partition(K, 0.9, count_mat, guide, sz, eleM, eleD, numP, repeat, maxiter, rank) ;
%end