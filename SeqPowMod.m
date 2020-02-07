function s = SeqPowMod(e, seq, P)

s=zeros(1, length(seq));
for i=1:length(seq)
    s(i)=powermod(e, seq(i), P);
end


end