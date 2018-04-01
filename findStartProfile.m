function TStart=findStartProfile(TP,P)

[maxP,imaxP]=max(P);

FSP=0.05*maxP; % Initial force value- find 5% from max force
iFSP=find(P(1:imaxP)<FSP); % find the indices before the max vel in which the velocity is bellow 5% max velocity    
if iFSP
    TStart=TP(iFSP(end)); % the last value in iFSP, is the index of the time in which force rises
else
    TStart=TP(1); % the last value in iFSP, is the index of the time in which force rises
end