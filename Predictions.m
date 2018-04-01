function Predictions

clc
clear 
close all

c = get(gcf,'DefaultAxesColorOrder');
ofs = get(gcf,'Position'); % original figure size
close

dF=0.4; % darkening factor
cnd=c(1,:)*1.2;
cnd_act=dF*cnd;
cd85=[247 181 106]/255; % light orange
cd85_act=dF*cd85;
cP=c(2,:)*1.167;
cV=c(5,:)*0.77;
cA=c(4,:)*0.9;
cVt=c(6,:)*0.4;

lw=8;
lwRep=6;
afs=19;
lfs=18;

sfs=1.3; % scale figure size
mfs=3; % move figure
fs=[ofs(1:2)/mfs ofs(3:4)*sfs];% figure size

yLim=[-1.5 5];

pi=0; % m
pf=0.1; % m

st=.005; % sample time
t=0:st:.55; % msec

delayT85=.085;
delayInd85=floor(delayT85/st);

state_model=0; % 0= from regression; 1= from taylor

% Pos
p=normcdf(linspace(0,1,length(t)),0.5,0.1)*(pf-pi)+pi;
% Vel
v=diff(p)./diff(t);
% Acc
a=diff(v)./diff(t(1:end-1));
% Delay vel
v_del85=[zeros(1,delayInd85),v(1:end-delayInd85)]; % for illustration - between 70 and 100 ms
% PertF
pertF_ND=v*6; % the perturbation force
pertF_D85=v_del85*6; % the perturbation force

t=t(1:end-2);
p=p(1:end-2);
v=v(1:end-1);
v_del85=v_del85(1:end-1);
pertF_ND=pertF_ND(1:end-1);
pertF_D85=pertF_D85(1:end-1);

% PV_ND- find the weights of pos and vel from the regression
b_PV_ND=regress(pertF_ND',[p' v']);
K_PV_ND=b_PV_ND(1)+4; % after playing with the values
B_PV_ND=b_PV_ND(2)-1.25;

% PV_ND- find the weights of pos and vel from the regression
b_PVt_D=regress(pertF_D85',[p' v_del85']);
K_PVt_D=b_PVt_D(1)+4; % after playing with the values
B_PVt_D=b_PVt_D(2)-1.25;

% PVA_D model- find the weights of pos, vel and acc from the regression
b_PVA_D=regress(pertF_D85',[p' v' a']);
K_PVA_D_reg=b_PVA_D(1)-5; % after playing with the values
B_PVA_D_reg=b_PVA_D(2)+1;
M_PVA_D_reg=b_PVA_D(3)+0.1;

if state_model
    K_PVA_D=K_PVA_D_tay;
    B_PVA_D=B_PVA_D_tay;
    M_PVA_D=M_PVA_D_tay;
else
    K_PVA_D=K_PVA_D_reg;
    B_PVA_D=B_PVA_D_reg;
    M_PVA_D=M_PVA_D_reg;
end

Act_ND=K_PV_ND*p+B_PV_ND*v;
Act_D_Time=K_PVt_D*p+B_PVt_D*v_del85;
Act_D_State=K_PVA_D*p+B_PVA_D*v+M_PVA_D*a;

t=t*1000;

% plot PertF vs time
figure('position',fs)
hold on
p4=plot(t,v*2,':','color',[1 1 1]*0.5,'linewidth',6);
p1=plot(t,pertF_ND,'color',cnd,'linewidth',lw);
p2=plot(t,pertF_D85,'color',cd85,'linewidth',lw);
line([0  t(end)] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
xlabel('Time [ms]')
ylabel('Perturbation Force [N]')
set(gca,'xtick',0:100:500,'xticklabel',-100:100:400,'fontsize',afs,'fontweight','b','box','off')
xlim([t(1) t(end)])
ylim([yLim(1) yLim(2)])
box off
l=legend([p4,p1,p2],' Velocity',' Non-delay',' Delay');
set(l,'fontweight','b','fontsize' ,lfs,'location','northwest');
legend boxoff

% plot PertF vs time
figure('position',fs)
hold on
p1=plot(t,pertF_ND,'color',cnd,'linewidth',lw);
p3=plot(t,K_PV_ND*p,':','color',cP,'linewidth',lw);
p4=plot(t,B_PV_ND*v,':','color',cV,'linewidth',lw);
p2=plot(t,Act_ND,'color',cnd_act,'linewidth',lwRep);
line([0  t(end)] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
xlabel('Time [ms]')
ylabel('Force [N]')
set(gca,'xtick',0:100:500,'xticklabel',-100:100:400,'fontsize',afs,'fontweight','b','box','off')
xlim([t(1) t(end)])
ylim([yLim(1) yLim(2)])
box off
l=legend([p1,p2,p3,p4],' Perturbation',' Representation',' Pos',' Vel');
set(l,'fontweight','b','fontsize' ,lfs,'location','northwest');
legend boxoff

% plot PertF vs time
figure('position',fs)
hold on
p1=plot(t,pertF_D85,'color',cd85,'linewidth',lw);
p3=plot(t,K_PVt_D*p,':','color',cP,'linewidth',lw);
p4=plot(t,B_PVt_D*v_del85,':','color',cVt,'linewidth',lw);
p2=plot(t,Act_D_Time,'color',cd85_act,'linewidth',lwRep);
line([0  t(end)] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
xlabel('Time [ms]')
ylabel('Force [N]')
set(gca,'xtick',0:100:500,'xticklabel',-100:100:400,'fontsize',afs,'fontweight','b','box','off')
xlim([t(1) t(end)])
ylim([yLim(1) yLim(2)])
box off
l=legend([p1,p2,p3,p4],' Perturbation',' Representation',' Pos',' Delayed Vel');
set(l,'fontweight','b','fontsize' ,lfs,'location','northwest');
legend boxoff

% plot PertF vs time
figure('position',fs)
hold on
p1=plot(t,pertF_D85,'color',cd85,'linewidth',lw);
p3=plot(t,K_PVA_D*p,':','color',cP,'linewidth',lw);
p4=plot(t,B_PVA_D*v,':','color',cV,'linewidth',lw);
p5=plot(t,M_PVA_D*a,':','color',cA,'linewidth',lw);
p2=plot(t,Act_D_State,'color',cd85_act,'linewidth',lwRep);
line([0  t(end)] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
xlabel('Time [ms]')
ylabel('Force [N]')
set(gca,'xtick',0:100:500,'xticklabel',-100:100:400,'fontsize',afs,'fontweight','b','box','off')
xlim([t(1) t(end)])
ylim([yLim(1) yLim(2)])
box off
l=legend([p1,p2,p3,p4,p5],' Perturbation',' Representation',' Pos',' Vel',' Acc');
set(l,'fontweight','b','fontsize' ,lfs,'location','northwest');
legend boxoff