%function run_partition(K, skewness, count_mat, guide, sz, eleM, eleD, numP, repeat, maxiter, rank)
% 1. Separate the dataset into K partitions
% 2. Securely aligning the element set
% 3. Run federated TF

%% partitioning
% get partitioning cutoff
if (skewness ==0) %balance. partition equal size
    cutoffs=1:floor(numP/K):(numP); %sz(2)= #patients, assuming patients 1:Ip
    cutoffs(K+1)=numP+1; %dummy for simplicity..
else %unbalanced data
    %Largest one take skenwness*numP, others share the rest. 
    cutoffs=[1,floor(numP*skewness):floor((numP-floor(numP*skewness))/(K-1)):(numP)];
    cutoffs(K+1)=numP+1;
end

% Save elements and update index for partitioned data (only for simulation.  Real data already have it)
eleMs=cell(1,K); %contains medi elements of kth partitioned data
eleDs=cell(1,K); %contains diag elements of kth partitoined data
idxs=cell(1,K);
%prevIdxsM=cell(1,K);
%prevIdxsD=cell(1,K);


for k=1:K
    %count_mat [count, patients, medi, diag, institution id]    
    idxs{1,k}=count_mat(:,2)>=cutoffs(k) & count_mat(:,2)<cutoffs(k+1);    
    %prevIdxsM{1,k}=count_mat(idxs{1,k}, 3);  
    %prevIdxsD{1,k}=count_mat(idxs{1,k}, 4);
    
    count_mat_k=count_mat(idxs{1,k}, :);
    idxM=unique(count_mat_k(:, 3))';
    idxD=unique(count_mat_k(:, 4))';
    %save elements
    eleMs{k}=double(eleM(1, idxM)); %medi
    eleDs{k}=double(eleD(1, idxD)); %diag
    %update index
    %newidxM=[idxM; 1:length(idxM)];
    %newidxD=[idxD; 1:length(idxD)];
    for i=1:length(idxM)
        count_mat_k(count_mat_k(:, 3)==idxM(i), 3)=i;
    end
    
    for i=1:length(idxD)
        count_mat_k(count_mat_k(:,4)==idxD(i), 4)=i;
    end
    
    count_mat(idxs{1,k}, :)=count_mat_k;
    
end
clearvars count_mat_k idxM idxD

%{
%Permutation approach
%add last column for institution id. 1:K
count_mat=[count_mat, zeros(size(count_mat,1),1)];
numP=sz(2);
% get partitioning cutoff
if (skewness ==0) %balance. partition equal size
    cutoffs=1:floor(numP/K):(numP); %sz(2)= #patients, assuming patients 1:Ip
    cutoffs(K+1)=numP+1; %dummy for simplicity..
else %unbalanced data
    %Largest one take skenwness*numP, others share the rest. 
    cutoffs=[1,floor(numP*skewness):floor((numP-floor(numP*skewness))/(K-1)):(numP)];
    cutoffs(K+1)=numP+1;
end

% permuting when patients at institutions have different medi and diag feature set
per_patients=randperm(sz(2));
%{
for p=1:sz(2)
    count_mat(count_mat(:,2)==p, 2)=per_patients(p);
end
count_mat(:,2)=count_mat(:,2)-sz(2);
clearvars per_patients 
%}
% partitioning
patients=cell(1, K);
for k=1:K
    patients{1,k}=per_patients(1,cutoffs(k):(cutoffs(k+1)-1));
    for ii=1:length(patients{1,k})
        count_mat(count_mat(:,2)==patients{1,k}(ii), end)=k;
    end
    %idx=count_mat(:,2)>=cutoffs(k) & count_mat(:,2)<cutoffs(k+1);
    %count_mat(idx,end)=k;
end

%FOR TESTING WHETHER the partition and merge is successful or not
count_mat_ori=count_mat; %before partitioned data are saved in count_mat_ori

% Save elements and update index for partitioned data (only for simulation.  Real data already have it)
eleMs=cell(1,K); %contains medi elements of kth partitioned data
eleDs=cell(1,K); %contains diag elements of kth partitoined data

for k=1:K
    %count_mat [count, patients, medi, diag, institution id]
    count_mat_k=count_mat(count_mat(:,5)==k, :);
    idxM=unique(count_mat_k(:, 3))';
    idxD=unique(count_mat_k(:, 4))';
    %save elements
    eleMs{k}=eleM(1, idxM); %medi
    eleDs{k}=eleD(1, idxD); %diag
    %update index
    %newidxM=[idxM; 1:length(idxM)];
    %newidxD=[idxD; 1:length(idxD)];
    for i=1:length(idxM)
        count_mat_k(count_mat_k(:,3)==idxM(i), 3)=i;
    end
    
    for i=1:length(idxD)
        count_mat_k(count_mat_k(:,4)==idxD(i), 4)=i;
    end
    
    count_mat(count_mat(:,5)==k, :)=count_mat_k;
    
end
%now count_mat contains new medi diag idx after partitioning 
clearvars count_mat_k idxM idxD eleM eleD
%}

%% Feature set alignment with sparse representation
sz_medi=4000; sz_diag=1000; % an arbitarilarly large number greater than the size of feature set.
P=randseed(1, 1,1, max([sz_medi, sz_diag])); % arbitary prime number larger than size of medi, diag
%find the updated idx
tAlign=tic;
[codeBooksM, codeBooksD, communication_align]=alignDimension(K, eleMs, eleDs, sz_medi, sz_diag, P);

newIdxsM=cell(1,K);
newIdxsD=cell(1,K);
for k=1:K
    %medi
    cur_CdBk_Srt=sortrows(codeBooksM{1, k},3); %sort accrd to prevIdx
    
    %prevIdx=prevIdxsM{1,k};
    prevIdx=count_mat(idxs{1,k}, 3); %medi
    newIdx=cur_CdBk_Srt(prevIdx,4);
    newIdxsM{1,k}=newIdx;
   % count_mat(count_mat(:, 5)==k, 3)=newIdx;
    
    %diag
    cur_CdBk_Srt=sortrows(codeBooksD{1, k},3); %sort accrd to prevIdx
    %prevIdx=prevIdxsD{1,k};
    prevIdx=count_mat(idxs{1,k}, 4); %diag
    newIdx=cur_CdBk_Srt(prevIdx,4);
    newIdxsD{1,k}=newIdx;
   % count_mat(count_mat(:, 5)==k, 4)=newIdx;
  
end

%clearvars prevIdxsM newIdxsD cur_CdBk_Srt prevIdx

%update new idx
for k=1:K
    count_mat(idxs{1,k}, 3)=newIdxsM{1,k};
    count_mat(idxs{1,k}, 4)=newIdxsD{1,k};
end
t_align=toc(tAlign);
clearvars newIdxsM newIdxsD prevIdxsM prevIdxsD


%% Sparse representation to sptensor
%Full data (not partitioned data)
%sptensor([coordinates], count, size)
Xs=cell(1,K);
%zeroIdxs=cell(3,K);
for k=1:K
    count_mat_k=count_mat(idxs{1,k},:);    
    Xs{k}=sptensor(count_mat_k(:, 2:4), count_mat_k(:,1), sz(2:end));
    
end


% Sparse tensor of full data
X=sptensor(count_mat(:,2:4), count_mat(:, 1), sz(2:end));


paraUpdateRatio=1.1;
lambda_q=1e-2;
omega=1;
gamma=1;


%% Run federated TF
results=cell(repeat, 6); %rmse, t, communication, t_iter, A, u
for r=1:repeat
    
    [results{r, 1},  results{r, 2}, results{r, 3}, results{r, 4}, results{r, 5}, results{r, 6}]=federatedTF(maxiter, Xs, X, sz(2:end), rank, K, guide, cutoffs,  omega, gamma, paraUpdateRatio, lambda_q);  

end


%% Save the results
%fileName=strcat('iteration_K',num2str(K),'skewness', num2str(skewness), '.mat');
%save(fileName);


%end
