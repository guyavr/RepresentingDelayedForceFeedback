function LogL=LogLikeLinearRegress(T,k,MSE)
% Calculates the Log Likelihood of a linear regression.
% According to http://www.le.ac.uk/users/dsgp1/COURSES/MATHSTAT/13mlreg.pdf
DOF=T-k; % Degrees of Freedom
LogL=-(T/2)*log(2*pi)-(T/2)*log(MSE)-0.5*DOF;