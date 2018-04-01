function ci = NonNanDataCI(data)
% This function takes the data matrix with NaN values, remove the Nan, and
% return 95% confidence interval for each column in data

nboot=1000;
bootfun=@(x)mean(x);
    
dataNoNan=data;
dataNoNan(isnan(data(:,1)),:)=[];
ci=bootci(nboot,bootfun,dataNoNan);

end

