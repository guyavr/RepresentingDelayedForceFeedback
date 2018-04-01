function [LogL, Npar]=parAICBIC(reg_b,reg_stats,T) 
% parAICBIC finds the parameters required for the aicbic function from the
% results of the linear regression (regress). T is the number of
% observation in the response vector of the regression
MSE=reg_stats(4); % Mean Squared Error
Npar=length(reg_b); % Number of parameters in the model.
LogL=LogLikeLinearRegress(T,Npar,MSE);