%Yejin Kim 2015 09 29
%Align dimension (Medication, Diagnosis) of seperated parties
%partitionedData=cell array of K parties seperated data
% partitoinMedication=matirx of each K parties' element in medication
% partitoinDiagnosis=matirx of each K parties' element in diagnosis

%K=number of institutions
%prtData=Dataset, sparse representation of tensor. Assume medi and diag are sorted
%increasing order. mode 1=patient, mode 2=medication, mode 3=diagnosis
%m2=total number of medication features (unknown, just guess) (or largest number of features, size of feature set)
%m3=total number of diagnosis features (unknown, just guess)
%prtM:feature set. always used as unique(prtM(k, :) and increasing sort). Not unique. Not increasing.for all
%institutions.
function [codeBookMs, codeBookDs, communication]= alignDimension(K, eleMs, eleDs, m2, m3, P)

%numRand=10;
%rmin=max(max(partitionsMedication))

aggrPolyM=zeros(K, m2+1); % for server
aggrPolyD=zeros(K, m3+1);
communication=0;

%% Distributed clients
%Make polynomials and multiplied with random polynomials
% Medications
%for k=1:K
for k=1:K
    
    polyM=ModPoly(leja(eleMs{1, k}),P);
    polyD=ModPoly(leja(eleDs{1, k}), P); %high degree polynomial is easy to give wrong value. Leja is sorting for better computation.
    
    polyM=[zeros(1,1+m2-length(polyM)),polyM];
    polyD=[zeros(1,1+m3-length(polyD)),polyD];
    
    %Send the polynomials to server
    aggrPolyD(k, :)=polyD; %% client
    aggrPolyM(k, :)=polyM; %% client
    
    %clearvars polyD polyM;
end
aggrPolyDBytes=whos('aggrPolyD');
aggrPolyMBytes=whos('aggrPolyM');
communication = communication + aggrPolyDBytes.bytes;
communication = communication + aggrPolyMBytes.bytes;

%% Server
% Aggregated polynomials
aggrPolyM(:, ~any(aggrPolyM, 1))=[];
aggrPolyD(:, ~any(aggrPolyD, 1))=[];
aggrPolyDBytes=whos('aggrPolyD');
aggrPolyMBytes=whos('aggrPolyM');
communication = communication + aggrPolyDBytes.bytes;
communication = communication + aggrPolyMBytes.bytes;

%{
aggrPolySmM=zeros([size(aggrPolyM), K]);
aggrPolySmD=zeros([size(aggrPolyD), K]);

for k=1:K
    aggrPolySmM(:, :, k)=mod(repmat(aggrPolyM(k, :), K,1) + aggrPolyM, P);
    aggrPolySmD(:, :, k)=mod(repmat(aggrPolyD(k, :), K,1) + aggrPolyD, P);
end

clearvars aggrPolyM aggrPolyD;
%}
% Send it back to all clients
% send aggrPolySumAtSrvMedi, aggrPolySumAtSrvDiag.

%% Distributed clients
% Find subset for all cases 2^(K-1)

aggrMmM=cell(K, 1); % For client to save membership of each elements.
aggrMmD=cell(K, 1);
aggrMmCntM=zeros(K, 2^K-1); % For server to aggregate the size
aggrMmCntD=zeros(K, 2^K-1); % For server to aggregate the size
%order of subset (1=included, 0=not included)
%001, 010 100 101 110 111
for k=1:K
%parfor k=1:K
  %[mmM, mmCntM]=findSubset(aggrPolySmM(:,:,k), eleMs{1, k}, k, K, P);
  [mmM, mmCntM]=findSubset(aggrPolyM, eleMs{1, k}, k, K, P);
  aggrMmM{k, :}= mmM;
  aggrMmCntM(k, :)=mmCntM; % send to server
  
  %[mmD, memCntD]=findSubset(aggrPolySmD(:,:,k), eleDs{1, k}, k, K, P);
  [mmD, memCntD]=findSubset(aggrPolyD, eleDs{1, k}, k, K, P);
  aggrMmD{k, :}= mmD;
  aggrMmCntD(k, :)=memCntD; % send to server
  
end
aggrMmMBytes=whos('aggrMmM');
aggrMmDBytes=whos('aggrMmD');
aggrMmCntMBytes=whos('aggrMmCntM');
aggrMmCntDBytes=whos('aggrMmCntD');
communication = communication + aggrMmMBytes.bytes;
communication = communication + aggrMmDBytes.bytes;
communication = communication + aggrMmCntMBytes.bytes;
communication = communication + aggrMmCntDBytes.bytes;

clearvars aggrPolySmM aggrPolySmD;

%% Server
% Send it back to all clients
aggrMmCntM=sum(aggrMmCntM);
% Duplicated count
for s=1:(2^K-1)
    aggrMmCntM(1, s)=aggrMmCntM(1, s)/sum(dec2bin(s)=='1');
end

aggrMmCntD=sum(aggrMmCntD);
for s=1:(2^K-1)
    aggrMmCntD(1, s)=aggrMmCntD(1, s)/sum(dec2bin(s)=='1');
end

%% Distributed clients
% Add zero rows to align
codeBookMs=cell(1,K);
codeBookDs=cell(1,K);

for k=1:K
%parfor k=1:K    
    %curCountM=count(count(:, 5)==k, 3);

    elements=eleMs{1,k}';
    codeBookM=[aggrMmM{k,:},  elements, (1:length(elements))', zeros(size(elements, 1), 1)]; %last column for updated idx
    codeBookM=sortrows(codeBookM); %aggrMmM=membership to each subset
    cumAggrMmCntM=cumsum(aggrMmCntM);
    cumAggrMmCntM=[0, cumAggrMmCntM]; %dummy 0
    eIdx=0;
    for s= unique(codeBookM(:,1))'
        sIdx=eIdx+1;
        eIdx=eIdx+aggrMmCntM(1,s);
        codeBookM(sIdx:eIdx, end)=((cumAggrMmCntM(1, s)+1):cumAggrMmCntM(1, s+1))';
        %newItemM((cumAggrMmCntM(1, s)+1):cumAggrMmCntM(1, s+1))=curItemM(sIdx:eIdx);
    end
    codeBookMs{1, k}=codeBookM;
    
    elements=eleDs{1, k}';
    codeBookD=[aggrMmD{k,:},  elements, (1:length(elements))', zeros(size(elements, 1), 1)]; %last column for updated idx
    codeBookD=sortrows(codeBookD); %aggrMmM=membership to each subset 
    cumAggrMmCntD=cumsum(aggrMmCntD);
    cumAggrMmCntD=[0, cumAggrMmCntD]; %dummy 0
    eIdx=0;
    for s= unique(codeBookD(:,1))'
        sIdx=eIdx+1;
        eIdx=eIdx+aggrMmCntD(1,s);
        codeBookD(sIdx:eIdx, end)=((cumAggrMmCntD(1, s)+1):cumAggrMmCntD(1, s+1))';
    end
    codeBookDs{1, k}=codeBookD;
    
    
    
    %{
    %Medication
    curData=prtData{k,:};
    %elements is a kinda code book? 
    elements=eleMs{1, k}';
    elements=[elements, (1:length(elements))'];
    elements=sortrows([aggrMmM{k,:},elements]); %group and sort elements
    
    curData=curData(:, elements(:,3)',:);%sort of non-zero Medications
    curItemM=elements(:,2)';
    
    [p, m, d]=size(curData);
    newCurData=zeros(p, sum(aggrMmCntM), d);
    newItemM=zeros(1, sum(aggrMmCntM));
    cumAggrMmCntM=cumsum(aggrMmCntM);
    cumAggrMmCntM=[0, cumAggrMmCntM]; %dummy 0
    eIdx=0;
    for s= unique(elements(:,1))'
        sIdx=eIdx+1;
        eIdx=eIdx+aggrMmCntM(1,s);
        newCurData(:, (cumAggrMmCntM(1, s)+1):cumAggrMmCntM(1, s+1),:)= curData(:, sIdx:eIdx, :);
        newItemM((cumAggrMmCntM(1, s)+1):cumAggrMmCntM(1, s+1))=curItemM(sIdx:eIdx);
    end
    
    
    curData=newCurData;
    
    %Diagnosis         
    elements=eleDs{1, k}';
    elements=[elements, (1:length(elements))'];
    elements=sortrows([aggrMmD{k,:},elements]); %group and sort elements
    
    curData=curData(:,:, elements(:,3)');%sort of non-zero Diagnosis
    curItemD=elements(:,2)';
    
    [p, m, d]=size(curData);
    newCurData=zeros(p, m, sum(aggrMmCntD));
    newItemD=zeros(1, sum(aggrMmCntD));
    cumAggrMmCntD=cumsum(aggrMmCntD);
    cumAggrMmCntD=[0, cumAggrMmCntD]; %dummy 0
    eIdx=0;
    for s= unique(elements(:,1))'
        sIdx=eIdx+1;
        eIdx=eIdx+aggrMmCntD(1,s);
        newCurData(:, :,(cumAggrMmCntD(1, s)+1):cumAggrMmCntD(1, s+1))= curData(:, :, sIdx:eIdx);
        newItemD((cumAggrMmCntD(1, s)+1):cumAggrMmCntD(1, s+1))=curItemD(sIdx:eIdx);
    end
    
    algnData{k,1}=newCurData;
    itemM(k, :)=newItemM;
    itemD(k, :)=newItemD;
    %}
end
aggrMmCntDBytes=whos('aggrMmCntD');
aggrMmCntMBytes=whos('aggrMmCntM');
communication = communication + aggrMmCntDBytes.bytes;
communication = communication + aggrMmCntMBytes.bytes;





end


