function [membership, mmCnt]=findSubset(evalMat, elements, k, K, P)
    %evalMat=aggrPolySm(:,:,k);
    %elements=unique(partitions(k, :));
    subsetsPair=zeros(K, length(elements));
    
    % find pairwise subset
    for kk= [1:k-1, k+1:K] % compare with kk
        coef=evalMat(kk,:); 
        subsetsPair(kk, (polyvalMod(coef, elements, P) == 0))=1;
        %{
        for e=1:length(elements)
            %{
            %geoSequence=elements(e).^[0:(length(coef)-1)];
            geoSequence=SeqPowMod(elements(e), 0:(length(coef)-1), P);
            geoSequence=fliplr(geoSequence);
            %geoSequence=mod(geoSequence, P);
            answers=mod(sum(mod(coef.*geoSequence, P)), P);
            %}
            answers=polyvalMod(coef, elements(e), P);
            if (answers==0)
                subsetsPair(kk, e)=1;
            end
        end
        %}
    end
    subsetsPair(k, :)=1;
    
    % find full subsets
    subsetsPair=subsetsPair';
    membership=bin2dec(num2str(subsetsPair)); %each elements' memberships
    %[count, sets]=hist(membership, unique(membership));
    [count, sets]=hist(membership, 1:(2^K-1));
    
    % Send the size of all subsets &  Aggregate the sizes at Server
    mmCnt=zeros(1, 2^K-1);
    mmCnt(sets)=count;
end