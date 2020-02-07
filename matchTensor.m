%Input: two factor matrices 1 by 3 cell
%A1:baseline, A2:ordered based on A1. 
%column order of A2 match to A1.
%output: sum of tensor after aligning phenotypes
function A2=matchTensor(A1, A2, rank)

%compute the cosine similarity of column of factor matrices at each mode
cost=zeros(rank, rank);
for r=1:rank
    for rr=1:rank
        cost_n=0;
        for n=2:3
            cost_n=cost_n+ A1{n}(:,r)'*A2{n}(:,rr)/(norm(A1{n}(:,r))*norm(A2{n}(:,rr)));
        end
        cost(r, rr)= cost_n;
    end
end

%similarity to dissimilarity
cost=-cost;
%find assignment using hungarian method
[match,]=hungarian(cost);
for n=1:3
    A2{n}=A2{n}(:,match);
end




end
