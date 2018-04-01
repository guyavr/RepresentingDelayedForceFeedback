function RepresentationAnalysis_D70SF_Generalization

clc
clear
close all

c = colormap(lines(10));
close
cS_subj=[70,130,180; 102,205,170; 64,224,208; 0,206,209; 32,178,170; 0,139,139; 135,206,235; 30,144,255]/255;
cF_subj=[238,130,238; 218,112,214; 147,112,219; 148,0,211; 139,0,139; 75,0,130; 123,104,238; 106,90,205]/255;

dF_rep=0.4; % darkening factor
dF_act=0.65; % darkening factor
scJb=0.5; % scale color for lightening jbfill

cV=c(5,:)*0.77;
cA=c(4,:)*0.9;
cVt=c(6,:)*0.4;

cS=c(6,:); % slow
cF=c(4,:)*1.2; % fast

s2ms=1000;

afs=19;
tfs=20;
lfs=18;

ofs = get(gcf,'Position'); % original figure size
close
sfs=1.3; % scale figure size
mfs=3; % move figure
fs=[ofs(1:2)/mfs ofs(3:4)*sfs];% figure size

lwRep=3;
lw=4;

xLim=[0 0.545]*s2ms;
yLim=[-0.6 2.8];
textloc=[0.35*s2ms 3.3]; % location of text box

ytick=[0:1:3];
xtick=[0:100:500];
xticklabel=[-100:100:400];

leg_act='Actual Force';

st=0.005;

[bwb,bwa]=butter(2,0.1); % build lowpass butterworth filter for the acceleration signal- 2nd order for 10Hz (10Hz/200Hz*2)

% Changing the data into force units by multiplying with the
% appropriate factors
VFactor=0.06; % N/cm
PFactor=0.21; % N*s/cm
AFactor=0.0082; % N*s^2/cm

dat=load('FC_BeFC_AllExpData_D70SF');
delVal=0.05;
delSamp=floor(delVal/st);
title='Group D70\_SF';

ReachVelMode=load('Reachmode.txt');
trVec_base=1:4; % force channel trials of only baseline
trVec_adapt=5:29; % force channel trials of only adaptation
trVec_wash=30:129; % force channel trials of only washout
Nt_wash=length(trVec_wash);

% plot forces of washout - separated for slow and fast according to reach
% vel mode
iWash=301:400;
RVM_wash=ReachVelMode(iWash);
% indices of slow and fast (according to instructions) from the beginning
% of washout
RVM_washSlow=find(RVM_wash==1);
RVM_washFast=find(RVM_wash==0);

% Split washout data to slow and fast (according to exp design and not actual velocities)
slow=dat.FC(:,trVec_adapt(end)+RVM_washSlow);
fast=dat.FC(:,trVec_adapt(end)+RVM_washFast);

trAna=1:5; % number of trials to analyze from the beginning of the washout for each velocity
Nt_ana=length(trAna);

Ns=size(dat.FC,1); % number of subjects

for velMode=1:2 % perform the analysis for each velocity: 1-slow, 2-fast
    if velMode==1
        FC=slow;
        cFC=cS;
        c_subj=cS_subj;
        name='Slow';
    else
        FC=fast;
        cFC=cF;
        c_subj=cF_subj;
        name='Fast';
    end
    
    c_act=cFC*dF_act;
    c_rep=cFC*dF_rep;
    
    % Find maximum force vector length to create an initialized vector
    TrialLen=nan(Ns,Nt_wash/2);
    
    for i=1:Nt_wash/2
        for s=1:Ns
            TrialLen(s,i)=length(FC(s,i).Force);
        end
    end
    
    % For only the analyzed trials
    TrialLenWash=reshape(TrialLen,[],1);
    
    Cut = prctile(TrialLenWash,10); % for 10% cut off
    lowCut=ceil(Cut(1));
    
    % Create matrices for TrialLen_FC and for TrialLen_BeFC that contain 1
    % for trials shorter than lowCut
    belowLowCutWash=zeros(Ns,Nt_wash/2);
    belowLowCutWash(TrialLen<lowCut)=1;
    
    for s=1:Ns
        
        for i=1:Nt_ana
            
            TrialLenAna=TrialLen(:,trAna);
            
            % For only the analyzed trials
            TrialLenAna_all=reshape(TrialLenAna,[],1);
            Cut = prctile(TrialLenAna_all,10); % for 10% cut off
            lowCut=ceil(Cut(1));
            
            % Create matrices for TrialLen_FC and for TrialLen_BeFC that contain 1
            % for trials shorter than lowCut
            belowLowCutWash=zeros(Ns,Nt_ana);
            belowLowCutWash(TrialLenAna<lowCut)=1;
            
            UniteTrLenAna(velMode)=lowCut-1; % It is -1 because of the acceleration vector
                        
        end
        
        
    end
    
    ActF.all=NaN(Ns*Nt_ana,UniteTrLenAna(velMode)); % take all the trAna trials of all subjects together
    PosY.all=NaN(Ns*Nt_ana,UniteTrLenAna(velMode)); % take all the trAna trials of all subjects together
    VelY.all=NaN(Ns*Nt_ana,UniteTrLenAna(velMode)); % take all the trAna trials of all subjects together
    VeltauY.all=NaN(Ns*Nt_ana,UniteTrLenAna(velMode)); % take all the trAna trials of all subjects together
    AccY.all=NaN(Ns*Nt_ana,UniteTrLenAna(velMode)); % take all the trAna trials of all subjects together
    
    ActF.mSubj=NaN(Ns,UniteTrLenAna(velMode)); % mean early generalization actual force for each subject
    PosY.mSubj=NaN(Ns,UniteTrLenAna(velMode)); % mean early generalization actual force for each subject
    VelY.mSubj=NaN(Ns,UniteTrLenAna(velMode)); % mean early generalization actual force for each subject
    VeltauY.mSubj=NaN(Ns,UniteTrLenAna(velMode)); % mean early generalization actual force for each subject
    AccY.mSubj=NaN(Ns,UniteTrLenAna(velMode)); % mean early generalization actual force for each subject
        
    for s=1:Ns
        for i=1:length(trAna)
            if ~belowLowCutWash(s,i)
                ActF.all((s-1)*Nt_ana+i,:)=FC(s,i).Force(1:UniteTrLenAna(velMode));
                PosY.all((s-1)*Nt_ana+i,:)=FC(s,i).PosY(1:UniteTrLenAna(velMode));
                VelY.all((s-1)*Nt_ana+i,:)=FC(s,i).VelY(1:UniteTrLenAna(velMode));
                AccY_org=diff(FC(s,i).VelY(1:UniteTrLenAna(velMode)+1))/st; % acceleration before filtering
                AccY.all((s-1)*Nt_ana+i,:)=filtfilt(bwb,bwa,AccY_org); % acceleration after filtering
                VeltauY.all((s-1)*Nt_ana+i,:)=[zeros(1,delSamp) FC(s,i).VelY(1:UniteTrLenAna(velMode)-delSamp)];
            end
        end
        ActF.mSubj(s,:)=nanmean(ActF.all((s-1)*Nt_ana+1:s*Nt_ana,:));
        PosY.mSubj(s,:)=nanmean(PosY.all((s-1)*Nt_ana+1:s*Nt_ana,:));
        VelY.mSubj(s,:)=nanmean(VelY.all((s-1)*Nt_ana+1:s*Nt_ana,:));
        VeltauY.mSubj(s,:)=nanmean(VeltauY.all((s-1)*Nt_ana+1:s*Nt_ana,:));
        AccY.mSubj(s,:)=nanmean(AccY.all((s-1)*Nt_ana+1:s*Nt_ana,:));
        
    end
    
    % Figures
    t=st*(1:UniteTrLenAna(velMode))*s2ms;
    
    ActF.m=nanmean(ActF.all); % mean of all subjects
    ActF.ci=NonNanDataCI(ActF.all);
    
    % Plot the ActForce
    figure('Position',fs)
    hold on
    p1=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    hold on
    h1=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    line([0  UniteTrLenAna(1)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    text(textloc(1),textloc(2),title,'fontsize' , tfs,'fontweight' , 'b')
    if velMode==1
        text(20,2.6,'Actual Force- Slow','fontsize' , lfs,'fontweight' , 'b')
    else
        text(20,2.6,'Actual Force- Fast','fontsize' , lfs,'fontweight' , 'b')
    end
    
    % Plot the ActForce for each subject
    figure('Position',fs)
    hold on
    for s=1:Ns
        p(s)=plot(t, ActF.mSubj(s,:),'-','color',c_subj(s,:), 'linewidth' ,lw); % data
    end
    line([0  UniteTrLenAna(1)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    if velMode==1
        text(20,2.6,'Actual Force- Slow','fontsize' , lfs,'fontweight' , 'b')
    else
        text(20,2.6,'Actual Force- Fast','fontsize' , lfs,'fontweight' , 'b')
    end
    
    if velMode==1
        ActF_Slow=ActF;
        PosY_Slow=PosY;
        VelY_Slow=VelY;
        VeltauY_Slow=VeltauY;
        AccY_Slow=AccY;
    else
        ActF_Fast=ActF;       
        PosY_Fast=PosY;
        VelY_Fast=VelY;
        VeltauY_Fast=VeltauY;
        AccY_Fast=AccY;
    end
end

tS=st*(1:UniteTrLenAna(1))*s2ms;
tF=st*(1:UniteTrLenAna(2))*s2ms;

% Test the generalization of the models: take the primitives gains from the
% slow movements during late adaptation and build the predicted actual forces for each
% participant according to the primitives in the fast movements
load('gs_D70SF_lateAdapt_PVA_PVVtau_VA_VVtau.mat');

% there is a decay in the actual forces from the alow movements of late
% adaptation (LA) and early washout (EW). Calculate the decay according to the ratio
% between the max actual force and max vel at LA.
load('ActF_VelY_factor_D70SF_LA.mat');

for velMode=1:2 % perform the analysis for each velocity: 1-slow, 2-fast
    if velMode==1
        FC=slow;
        cFC=cS;
        t=tS;
        name='Slow';
        
        ActF=ActF_Slow;        
        PosY=PosY_Slow;
        VelY=VelY_Slow;
        VeltauY=VeltauY_Slow;
        AccY=AccY_Slow;
    else
        FC=fast;
        cFC=cF;
        t=tF;
        name='Fast';
        
        ActF=ActF_Fast;        
        PosY=PosY_Fast;
        VelY=VelY_Fast;
        VeltauY=VeltauY_Fast;
        AccY=AccY_Fast;
    end
    
    c_act=cFC*dF_act;
    c_rep=cFC*dF_rep;
    
    maxVelY=max(nanmean(VelY.all));
    % calculate the predicted actual force for slow/fast movements according to
    % their ActF/VelY ratio from LA
    maxActF_predicted=maxVelY*ActF_VelY_factor_D70SF_LA;
    maxActF_actual=max(ActF.m);
    % The force decay is the ratio between the actual and predicted decay.
    LA2EW_forceDecay=maxActF_actual/maxActF_predicted;
    
    % VA
    GenF.VA.BV.all=NaN(Ns,length(t)); % mean for each subject
    GenF.VA.MA.all=NaN(Ns,length(t)); % mean for each subject
    
    for s=1:Ns
        GenF.VA.BV.all(s,:)=VFactor*gs_VA_lateAdapt(s,1)*VelY.mSubj(s,:)*LA2EW_forceDecay;
        GenF.VA.MA.all(s,:)=AFactor*gs_VA_lateAdapt(s,2)*AccY.mSubj(s,:)*LA2EW_forceDecay;
        GenF.VA.Tot.all=GenF.VA.BV.all + GenF.VA.MA.all;
    end
    
    GenF.VA.BV.m=nanmean(GenF.VA.BV.all); % mean of all subjects
    GenF.VA.BV.ci=NonNanDataCI(GenF.VA.BV.all);
    GenF.VA.MA.m=nanmean(GenF.VA.MA.all); % mean of all subjects
    GenF.VA.MA.ci=NonNanDataCI(GenF.VA.MA.all);
    GenF.VA.Tot.m=nanmean(GenF.VA.Tot.all); % mean of all subjects
    GenF.VA.Tot.ci=NonNanDataCI(GenF.VA.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p4=plot(t, GenF.VA.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    p5=plot(t, GenF.VA.MA.m,':', 'color' ,cA , 'linewidth' ,lw);
    p1=plot(t, GenF.VA.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,GenF.VA.Tot.ci(2,:),GenF.VA.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,GenF.VA.BV.ci(2,:),GenF.VA.BV.ci(1,:),1-((1-cV)*scJb));
    jbfill(t,GenF.VA.MA.ci(2,:),GenF.VA.MA.ci(1,:),1-((1-cA)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h4=plot(t, GenF.VA.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    h5=plot(t, GenF.VA.MA.m,':', 'color' ,cA , 'linewidth' ,lw);
    h1=plot(t, GenF.VA.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  t(end)] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h4,h5], leg_act, 'Predicted Generalization' ,'Vel','Acc');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    text(380,2.57,name,'fontsize' , tfs,'fontweight' , 'b')
    
    % VVtau
    GenF.VVtau.BV.all=NaN(Ns,length(t)); % mean for each subject
    GenF.VVtau.BtVtau.all=NaN(Ns,length(t)); % mean for each subject
    
    for s=1:Ns
        GenF.VVtau.BV.all(s,:)=VFactor*gs_VVtau_lateAdapt(s,1)*VelY.mSubj(s,:)*LA2EW_forceDecay;
        GenF.VVtau.BtVtau.all(s,:)=VFactor*gs_VVtau_lateAdapt(s,2)*VeltauY.mSubj(s,:)*LA2EW_forceDecay;
        GenF.VVtau.Tot.all=GenF.VVtau.BV.all + GenF.VVtau.BtVtau.all;
    end
    
    GenF.VVtau.BV.m=nanmean(GenF.VVtau.BV.all); % mean of all subjects
    GenF.VVtau.BV.ci=NonNanDataCI(GenF.VVtau.BV.all);
    GenF.VVtau.BtVtau.m=nanmean(GenF.VVtau.BtVtau.all); % mean of all subjects
    GenF.VVtau.BtVtau.ci=NonNanDataCI(GenF.VVtau.BtVtau.all);
    GenF.VVtau.Tot.m=nanmean(GenF.VVtau.Tot.all); % mean of all subjects
    GenF.VVtau.Tot.ci=NonNanDataCI(GenF.VVtau.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p4=plot(t, GenF.VVtau.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    p5=plot(t, GenF.VVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    p1=plot(t, GenF.VVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,GenF.VVtau.Tot.ci(2,:),GenF.VVtau.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,GenF.VVtau.BV.ci(2,:),GenF.VVtau.BV.ci(1,:),1-((1-cV)*scJb));
    jbfill(t,GenF.VVtau.BtVtau.ci(2,:),GenF.VVtau.BtVtau.ci(1,:),1-((1-cVt)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h4=plot(t, GenF.VVtau.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    h5=plot(t, GenF.VVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    h1=plot(t, GenF.VVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  t(end)] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h4,h5], leg_act, 'Predicted Generalization' ,'Vel','Delayed Vel');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    text(380,2.57,name,'fontsize' , tfs,'fontweight' , 'b')
end