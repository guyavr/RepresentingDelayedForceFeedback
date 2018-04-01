function Predictions_Generalization_slow2fast
% Generalization to different velocities: after adapting to to delayed
% velocity-dependent force field with movement duration of .55 and
% constructing the primitives, what would be the force during sudden trials
% with force channels, in which the the movement duration is larger
% (slower) or smaller (faster).

clc
clear
close all

c = get(gcf,'DefaultAxesColorOrder');
ofs = get(gcf,'Position'); % original figure size
close

cP=c(2,:)*1.167;
cV=c(5,:)*0.77;
cA=c(4,:)*0.8;
cVt=c(6,:)*0.4;

cS=c(6,:); % slow
cF=c(4,:)*1.2; % fast

lw=8;
lwRep=6;
afs=19;
lfs=18;

sfs=1.3; % scale figure size
mfs=3; % move figure
fs=[ofs(1:2)/mfs ofs(3:4)*sfs];% figure size

yLim=[-1.2 2.2];
yTick=-1:3;

pi=0; % cm
pf=10; % cm

st=.005; % sample time

delayT=.07;
delayInd=floor(delayT/st);

tVel=[.55, .3]; % t for slow movements, original movement velocity (first exp), and fast

% Changing the data into force units by multiplying with the
% appropriate factors
PFactor=0.32; % N/cm
VFactor=0.06; % N*s/cm
AFactor=0.0068; % N*s^2/cm

% the weights of the primitives for each representation model (PVA, PVVt)
% following adaptation to delayed force in .55 movement duration (taken
% from the experimental results, according to the mean of each).
load('g_D70_lateAdapt_PVA_PVVtau')

% multiply by a factor that would decrease the forces - because the
% movement during adaptation are slower (upper duration is smaller than the
% D70 group). Also, there is a decay in force throughout washout.
K_PVA=g_PVA_lateAdapt(1)*0.3;
B_PVA=g_PVA_lateAdapt(2)*0.3;
M_PVA=g_PVA_lateAdapt(3)*0.3;

K_PVVt=g_PVVtau_lateAdapt(1)*0.3;
B_PVVt=g_PVVtau_lateAdapt(2)*0.3;
Bt_PVVt=g_PVVtau_lateAdapt(3)*0.3;

for v=1:length(tVel)
    
    if v==1
        cVel=cS;
    else
        cVel=cF;
    end
    t=0:st:tVel(v); % sec
    % Pos
    p=normcdf(linspace(0,1,length(t)),0.5,0.1)*(pf-pi)+pi;
    % Vel
    v=diff(p)./diff(t);
    % Acc
    a=diff(v)./diff(t(1:end-1));
    % Delay vel
    v_del=[zeros(1,delayInd),v(1:end-delayInd)];
    
    t=t(1:end-2)*1000;
    p=p(1:end-2);
    v=v(1:end-1);
    v_del=v_del(1:end-1);
    
    Act_D_PVA=PFactor*K_PVA*p+VFactor*B_PVA*v+AFactor*M_PVA*a;
    Act_D_PVVt=PFactor*K_PVVt*p+VFactor*B_PVVt*v+VFactor*Bt_PVVt*v_del;
    
    % plot PVA vs time
    figure('position',fs)
    hold on
    p3=plot(t,PFactor*K_PVA*p,':','color',cP,'linewidth',lw);
    p4=plot(t,VFactor*B_PVA*v,':','color',cV,'linewidth',lw);
    p5=plot(t,AFactor*M_PVA*a,':','color',cA,'linewidth',lw);
    p2=plot(t,Act_D_PVA,'color',cVel,'linewidth',lwRep);
    line([0  tVel(1)]*1000 , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    xlabel('Time [ms]')
    ylabel('Force [N]')
    set(gca,'fontsize',afs,'xtick',0:100:500,'xticklabel',-100:100:400,'ytick',yTick,'fontweight','b','box','off')
    xlim([t(1) tVel(1)]*1000)
    ylim([yLim(1) yLim(2)])
    box off
    l=legend([p2,p3,p4,p5],' Predicted Actual Force',' Pos',' Vel',' Acc');
    set(l,'fontweight','b','fontsize' ,lfs,'location','northeast');
    legend boxoff
    
    % plot PVVt vs time
    figure('position',fs)
    hold on
    p3=plot(t,PFactor*K_PVVt*p,':','color',cP,'linewidth',lw);
    p4=plot(t,VFactor*B_PVVt*v,':','color',cV,'linewidth',lw);
    p5=plot(t,VFactor*Bt_PVVt*v_del,':','color',cVt,'linewidth',lw);
    p2=plot(t,Act_D_PVVt,'color',cVel,'linewidth',lwRep);
    line([0  tVel(1)]*1000 , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    xlabel('Time [ms]')
    ylabel('Force [N]')
    set(gca,'fontsize',afs,'xtick',0:100:500,'xticklabel',-100:100:400,'ytick',yTick,'fontweight','b','box','off')
    xlim([t(1) tVel(1)]*1000)
    ylim([yLim(1) yLim(2)])
    box off
    l=legend([p2,p3,p4,p5],' Predicted Actual Force',' Pos',' Vel',' Delayed Vel');
    set(l,'fontweight','b','fontsize' ,lfs,'location','northeast');
    legend boxoff
    
end