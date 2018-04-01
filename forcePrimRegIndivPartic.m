function varargout = forcePrimRegIndivPartic(f,prim,primdum,Ns)
% The function perform the regression between the actual force f and a
% combination of primitives. prim is the matrix of primitives with the
% dummy variables matrix. the function returns the regression fit R2,
% Log Likelihood (LogL), number of parameters in the regression model
% (Npar) and the primitives weights. the number of output weights is
% according to the number of primitives in the model.
nout = size(prim,2);
primprimdum=[prim primdum];

b=regress(f,primprimdum);
varargout=cell(1,nout);
for k = 1:nout
    varargout{k}=nan(Ns,1);
    ind_bS=nout+(Ns-1)*(k-1)+(1:(Ns-1));
    % weights
    varargout{k}(1:Ns-1) = b(k)+b(ind_bS);
    varargout{k}(Ns)= b(k)- sum(b(ind_bS));
end


