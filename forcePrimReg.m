function [R2,LogL,Npar,varargout] = forcePrimReg(f,prim,primdum)
% The function perform the regression between the actual force f and a
% combination of primitives. prim is the matrix of primitives with the
% dummy variables matrix. the function returns the regression fit R2,
% Log Likelihood (LogL), number of parameters in the regression model
% (Npar) and the primitives weights. the number of output weights is
% according to the number of primitives in the model.
nout = 2*size(prim,2);
primprimdum=[prim primdum];
T=length(f);
[b, bint, r, rint, stats]=regress(f,primprimdum);

R2=stats(1);
% calculating Log-Likelihood
[LogL, Npar]=parAICBIC(b,stats,T);
for k = 1:nout/2
    % weights
    varargout{2*(k-1)+1} = b(k);
    varargout{2*(k-1)+2} = bint(k,:);
end


