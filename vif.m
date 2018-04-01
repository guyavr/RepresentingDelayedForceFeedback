function V = vif(X)
% Variance Inflation Factor

nPred=size(X,2);
nObs=size(X,1);
V=nan(1,nPred);
for j=1:nPred
    jPred=setdiff(1:nPred,j);
    [~, ~, ~, ~, stats]=regress(X(:,j),[ones(nObs,1) X(:,jPred)]);
    R2=stats(1);
    V(j)=1/(1-R2);
end

end

