function prim = primDummy_differentTrLen( primConcat,primLen )
% Creates a matrix with a primitive after multiplication with dummy
% variables
% primConcat is a primitive concatenation for multiple trials and/or
% subjects.
% primLen is a vector, with a size of 1 X number of subjects, contaiing the
% length of each concatenated primitive for each subject.

sNum=length(primLen); % number of subjects being analyzed
prim=zeros(sum(primLen(1:end-1)),sNum-1);
for s=1:sNum
    indToFill=sum(primLen(1:s-1))+1:sum(primLen(1:s));
    if s<sNum
        prim(indToFill,s)=primConcat(indToFill);
    else
        lastSubj=repmat(primConcat(indToFill),1,sNum-1);
    end
end
prim=[prim; -1*lastSubj];
        
end

