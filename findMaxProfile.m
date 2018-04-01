function [maxP,tmaxP]=findMaxProfile(TP,P)
% TP is the Time vector of the Profile; P is the Profile
[maxP,imaxP]=max(P);
tmaxP=TP(imaxP);  
