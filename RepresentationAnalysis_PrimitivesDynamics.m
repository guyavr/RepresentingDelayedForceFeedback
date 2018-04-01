function RepresentationAnalysis_PrimitivesDynamics

clc
clear
close all

c = colormap(lines(10));
close
cnd_pert=c(1,:)*1.2;
cd70_pert=c(3,:)*1.05;
cd100_pert=c(7,:)*1.05;
scJb=0.5; % scale color for lightening jbfill

cP=c(2,:)*1.167;
cV=c(5,:)*0.77;
cVt=c(6,:)*0.4;
cm=linspace(0,1,25)';

afs=19;
tfs=20;
lfs=18;

ofs = get(gcf,'Position'); % original figure size
close
sfs=1.3; % scale figure size
mfs=3; % move figure
fstimedyn=[ofs(1:2)/mfs ofs(3)*1.2*sfs ofs(4)*sfs];% figure size
fs=[ofs(1:2)/mfs ofs(3)*sfs ofs(4)*sfs];% figure size
fsg=[ofs(1:2)/mfs ofs(3)*0.7*sfs ofs(4)*sfs]; % for gain
    
lw=2;
ms=9; % marker size

textloc=[6 0.9]; % location of text box

st=0.005;

[bwb,bwa]=butter(2,0.1); % build lowpass butterworth filter for the acceleration signal- 2nd order for 10Hz (10Hz/200Hz*2)

trVec_adapt=6:30; % force channel trials of only adaptation
Nt_adapt=length(trVec_adapt);

% Changing the data into force units by multiplying with the
% appropriate factors
VFactor=0.06; % N/cm

vifThr=10;

VIF_PV=cell(3,45);
VIF_VVt=cell(3,45);
VIF_PVVt=cell(3,45);

VIF_PV_K=nan(3,45);
VIF_PV_B=nan(3,45);

VIF_VVt_B=nan(3,45);
VIF_VVt_Bt=nan(3,45);

VIF_PVVt_K=nan(3,45);
VIF_PVVt_B=nan(3,45);
VIF_PVVt_Bt=nan(3,45);

VIF_SubjTr_PV_K=cell(3,1);
VIF_SubjTr_PV_B=cell(3,1);

VIF_SubjTr_PVVt_K=cell(3,1);
VIF_SubjTr_PVVt_B=cell(3,1);
VIF_SubjTr_PVVt_Bt=cell(3,1);

islargeVIF_PV=cell(3,1);
islargeVIF_PVVt=cell(3,1);

K_PV_trialSubj=cell(3,1);
B_PV_trialSubj=cell(3,1);

K_PVVt_trialSubj=cell(3,1);
B_PVVt_trialSubj=cell(3,1);
Bt_PVVt_trialSubj=cell(3,1);

for g=1:3 % 1- ND, 2- D70, 3- D100
    
    if g==1
        dat=load('FC_BeFC_AllExpData_ND');
        c_pert=cnd_pert;
        delVal=0;
        title='Group ND';
        PFactor=0.38; % N*s/cm
        AFactor=0.0056; % N*s^2/cm
    elseif g==2
        dat=load('FC_BeFC_AllExpData_D70');
        c_pert=cd70_pert;
        delVal=0.07;
        title='Group D70';
        PFactor=0.32; % N*s/cm
        AFactor=0.0068; % N*s^2/cm
    else
        dat=load('FC_BeFC_AllExpData_D100');
        c_pert=cd100_pert;
        delVal=0.1;
        title='Group D100';
        PFactor=0.26; % N*s/cm
        AFactor=0.0073; % N*s^2/cm
    end
    
    Nt=size(dat.FC,2);
    delSamp=floor(delVal/st);
    
    Ns=size(dat.FC,1); % number of subjects
    
    % Arrange actual force data according to trial- for each trial an array of
    % Ns (number of subjects) force profiles (in raws) - only for adaptation
    
    % Find maximum force vector length to create an initialized vector
    TrialLen_FC=nan(Ns,Nt);
    TrialLen_BeFC=nan(Ns,Nt);
    
    for i=1:Nt
        
        for s=1:Ns
            TrialLen_FC(s,i)=length(dat.FC(s,i).Force);
            TrialLen_BeFC(s,i)=length(dat.BeFC(s,i).Force);
        end
        
    end
    
    % Analyse the statistics of the duration of all FC and BeFC trials (for
    % all subjects)
    TrialLen_All=[reshape(TrialLen_FC,numel(TrialLen_FC),1); reshape(TrialLen_BeFC,numel(TrialLen_BeFC),1)];
    
    Cut = prctile(TrialLen_All,10); % for 10% cut off 
    lowCut=ceil(Cut(1));
    
    % Create matrices for TrialLen_FC and for TrialLen_BeFC that contain 1
    % for trials shorter than lowCut
    belowLowCut_FC=zeros(Ns,Nt);
    belowLowCut_BeFC=zeros(Ns,Nt);
    belowLowCut_FC(TrialLen_FC<lowCut)=1;
    belowLowCut_BeFC(TrialLen_BeFC<lowCut)=1;
    
    indLowCut=union(find(belowLowCut_FC),find(belowLowCut_BeFC));
    belowLowCut=zeros(Ns,Nt);
    belowLowCut(indLowCut)=1;
    
    TrBelowLowCut(g)=sum(sum(belowLowCut)); % trial pairs in each group that were below lowcut
    UniteTrLen(g)=lowCut-1; % It is -1 because of the acceleration vector
    
    maxVel=nan(Ns,Nt); % find max velocity for VFactor
    maxAcc=nan(Ns,Nt); % find max velocity for VFactor
    
    VIF_SubjTr_PV_K{g}=nan(Ns,Nt);
    VIF_SubjTr_PV_B{g}=nan(Ns,Nt);
    
    VIF_SubjTr_PVVt_K{g}=nan(Ns,Nt);
    VIF_SubjTr_PVVt_B{g}=nan(Ns,Nt);
    VIF_SubjTr_PVVt_Bt{g}=nan(Ns,Nt);
    
    islargeVIF_PV{g}=zeros(Ns,Nt);
    islargeVIF_PVVt{g}=zeros(Ns,Nt);
    
    K_PV_trialSubj{g}=nan(Ns,Nt);
    B_PV_trialSubj{g}=nan(Ns,Nt);
    
    K_PVVt_trialSubj{g}=nan(Ns,Nt);
    B_PVVt_trialSubj{g}=nan(Ns,Nt);
    Bt_PVVt_trialSubj{g}=nan(Ns,Nt);

    for i=1:Nt
        % To find a primitives weights of all subjects in one trial- concatenate the
        % forces and primitives of all subjects;
        Force_FC_trial=[];
        Force_BeFC_trial=[];
        PosY_trial=[];
        VelY_trial=[];
        AccY_trial=[];
        VelYtau_trial=[];
        
        trLenAnalyzeTrial=[];
        for s=1:Ns
            TrLen=min(TrialLen_FC(s,i),TrialLen_BeFC(s,i))-1;
            
            Force_FC=dat.FC(s,i).Force(1:TrLen)';
            
            PosY=dat.BeFC(s,i).PosY(1:TrLen)';
            VelY=dat.BeFC(s,i).VelY(1:TrLen)';
            AccY_org=diff(dat.BeFC(s,i).VelY(1:TrLen+1)')/st; % acceleration before filtering
            AccY=filtfilt(bwb,bwa,AccY_org); % acceleration after filtering
            VeltauY=[zeros(1,delSamp) dat.BeFC(s,i).VelY(1:TrLen-delSamp)]';
            
            if VelY
                % Calculate for only BeFC trials
                maxVel(s,i)=max(VelY);
                maxAcc(s,i)=max(abs(AccY));
                
                X_SubjTr_PV=[PosY VelY];
                VIF_SubjTr_PV=vif(X_SubjTr_PV);
                VIF_SubjTr_PV_K{g}(s,i)=VIF_SubjTr_PV(1);
                VIF_SubjTr_PV_B{g}(s,i)=VIF_SubjTr_PV(2);
                
                X_SubjTr_PVVt=[PosY VelY VeltauY];
                VIF_SubjTr_PVVt=vif(X_SubjTr_PVVt);
                VIF_SubjTr_PVVt_K{g}(s,i)=VIF_SubjTr_PVVt(1);
                VIF_SubjTr_PVVt_B{g}(s,i)=VIF_SubjTr_PVVt(2);
                VIF_SubjTr_PVVt_Bt{g}(s,i)=VIF_SubjTr_PVVt(3);
                
            end
            
            if VIF_SubjTr_PV_K{g}(s,i)> vifThr || VIF_SubjTr_PV_B{g}(s,i)> vifThr
                islargeVIF_PV{g}(s,i)=1;
            end
            
            if VIF_SubjTr_PVVt_K{g}(s,i)> vifThr || VIF_SubjTr_PVVt_B{g}(s,i)> vifThr || VIF_SubjTr_PVVt_Bt{g}(s,i)> vifThr
                islargeVIF_PVVt{g}(s,i)=1;
            end
            
            if g==1
                islargeVIF=islargeVIF_PV{g};
            else
                islargeVIF=islargeVIF_PVVt{g};
            end
            
            % for calculating the primitives's weight for each trial
            if TrLen>=UniteTrLen(g) && ~islargeVIF(s,i) % keep only movements that are not too fast
                Force_FC_trial=[Force_FC_trial;Force_FC];
                PosY_trial=[PosY_trial;PosY];
                VelY_trial=[VelY_trial;VelY];
                AccY_trial=[AccY_trial;AccY];
                VelYtau_trial=[VelYtau_trial;VeltauY];
                
                % To create the primitives matrix- should be with dummy variables
                % according to the subjects analyzed in each trial
                trLenAnalyzeTrial=[trLenAnalyzeTrial,TrLen]; % subjects to analyze is length(trLenAnalyzeTrial)
                
                % primitive gains for each subject in each trial - for
                % statistical analysis between the velocity and delayed
                % velocity gains between the groups
                % pos (P) and vel (V)
                prim=[PosY*PFactor VelY*VFactor];
                [b, ~, ~, ~, ~]=regress(Force_FC,prim);
                K_PV_trialSubj{g}(s,i)=b(1);
                B_PV_trialSubj{g}(s,i)=b(2);
                
                prim=[PosY*PFactor VelY*VFactor VeltauY*VFactor];
                [b, ~, ~, ~, ~]=regress(Force_FC,prim);
                K_PVVt_trialSubj{g}(s,i)=b(1);
                B_PVVt_trialSubj{g}(s,i)=b(2);
                Bt_PVVt_trialSubj{g}(s,i)=b(3);

            end
            
        end
        
        FAct=Force_FC_trial;
        P=PosY_trial*PFactor;
        V=VelY_trial*VFactor;
        A=AccY_trial*AFactor;
        Vtau=VelYtau_trial*VFactor;
        
        Pdum = primDummy_differentTrLen(P,trLenAnalyzeTrial);
        Vdum = primDummy_differentTrLen(V,trLenAnalyzeTrial);
        Vtaudum = primDummy_differentTrLen(Vtau,trLenAnalyzeTrial);
        Adum = primDummy_differentTrLen(A,trLenAnalyzeTrial);
        
        % pos (P) and vel (V)
        prim=[P V]; primdum=[Pdum Vdum];
        [R2_PV(g,i),LogL_PV,Npar_PV,K_PV(g,i),Kint_PV(g,i,:),B_PV(g,i),Bint_PV(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        X_PV=[prim primdum];
        VIF_PV{g,i}=vif(X_PV);
        VIF_PV_K(g,i)=VIF_PV{g,i}(1);
        VIF_PV_B(g,i)=VIF_PV{g,i}(2);
        
        % vel (V) and delayed vel (Vtau)
        prim=[V Vtau]; primdum=[Vdum Vtaudum];
        [R2_VVtau(g,i),LogL_VVtau,Npar_VVtau,B_VVtau(g,i),Bint_VVtau(g,i,:),Bt_VVtau(g,i),Btint_VVtau(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        X_VVt=[prim primdum];
        VIF_VVt{g,i}=vif(X_VVt);
        VIF_VVt_B(g,i)=VIF_VVt{g,i}(1);
        VIF_VVt_Bt(g,i)=VIF_VVt{g,i}(2);
        
        % pos (P), vel (V) and delayed vel (Vtau)
        prim=[P V Vtau]; primdum=[Pdum Vdum Vtaudum];
        [R2_PVVtau(g,i),LogL_PVVtau,Npar_PVVtau,K_PVVtau(g,i),Kint_PVVtau(g,i,:),B_PVVtau(g,i),Bint_PVVtau(g,i,:),Bt_PVVtau(g,i),Btint_PVVtau(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        X_PVVt=[prim primdum];
        VIF_PVVt{g,i}=vif(X_PVVt);
        VIF_PVVt_K(g,i)=VIF_PVVt{g,i}(1);
        VIF_PVVt_B(g,i)=VIF_PVVt{g,i}(2);
        VIF_PVVt_Bt(g,i)=VIF_PVVt{g,i}(3);
        
    end
    
    % mean maxVel and maxAcc of all BeFC trials
    mmMaxVel(g)=nanmean(nanmean(maxVel));
    mmMaxAcc(g)=nanmean(nanmean(maxAcc));
        
    % figures
    if g==1
        figure('Position',fstimedyn)
        hold on
        line([5.5 5.5],[-0.2 1], 'Color' , [0.6  0.6  0.6],'linestyle' , '--'  , 'linewidth' , 2);
        line([30.5 30.5],[-0.2 1], 'Color' , [0.6  0.6  0.6],'linestyle' , '--'  , 'linewidth' , 2);
        p1=plot(1:Nt,K_PV(g,:),'-o','color',cP,'markerfacecolor',cP,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        p2=plot(1:Nt,B_PV(g,:),'-^','color',cV,'markerfacecolor',cV,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        jbfill(1:Nt,Kint_PV(g,:,2),Kint_PV(g,:,1),1-((1-cP)*scJb));
        jbfill(1:Nt,Bint_PV(g,:,2),Bint_PV(g,:,1),1-((1-cV)*scJb));
        hold on
        h1=plot(1:Nt,K_PV(g,:),'-o','color',cP,'markerfacecolor',cP,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        h2=plot(1:Nt,B_PV(g,:),'-^','color',cV,'markerfacecolor',cV,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        line([1  Nt] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
        ylim([-0.2 1]);
        xlim([1  Nt]);
        text(textloc(1),textloc(2),title,'fontsize' , tfs,'fontweight' , 'b')
        xlabel('Force Channel Trial Number','fontsize' , 20, 'fontweight' , 'b');
        ylabel('Gain','fontsize' , 20, 'fontweight' , 'b');
        l=legend([h1,h2], 'Pos', 'Vel');
        set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northeast');
        legend boxoff
        set(gca,'box','off', 'FontSize', afs,'fontweight','b','ytick',0:0.2:1)
        
        figure('Position',fsg)
        hold on
        for i=1:Nt_adapt
            p(i)=plot(K_PV(g,trVec_adapt(i)),B_PV(g,trVec_adapt(i)),'.','markersize',30,'color',c_pert*cm(i),'linewidth',lw);
        end
        text(15,0.9,title,'fontsize' , tfs,'fontweight' , 'b')
        xlabel('Position Gain','fontsize' , 20, 'fontweight' , 'b');
        ylabel('Velocity Gain','fontsize' , 20, 'fontweight' , 'b');
        set(gca,'xtick',0:0.2:0.2,'box','off', 'FontSize', afs,'fontweight','b')
        axis equal
        ylim([-0.1 1]);
        xlim([-0.1 0.3]);
        trColMap=cm*c_pert;
        colormap(trColMap);
        hnd=colorbar;
        ylabel(hnd,'Force Channel Trial Number','fontsize',20);
        set(hnd,'box','off','ticks',[0.2:0.2:1],'yticklabel',5:5:25);
        
    else

        % PVVtau
        figure('Position',fstimedyn)
        hold on
        line([5.5 5.5],[-0.2 1], 'Color' , [0.6  0.6  0.6],'linestyle' , '--'  , 'linewidth' , 2);
        line([30.5 30.5],[-0.2 1], 'Color' , [0.6  0.6  0.6],'linestyle' , '--'  , 'linewidth' , 2);
        p1=plot(1:Nt,K_PVVtau(g,:),'-o','color',cP,'markerfacecolor',cP,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        p2=plot(1:Nt,B_PVVtau(g,:),'-^','color',cV,'markerfacecolor',cV,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        p3=plot(1:Nt,Bt_PVVtau(g,:),'-s','color',cVt,'markerfacecolor',cVt,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        jbfill(1:Nt,Kint_PVVtau(g,:,2),Kint_PVVtau(g,:,1),1-((1-cP)*scJb));
        jbfill(1:Nt,Bint_PVVtau(g,:,2),Bint_PVVtau(g,:,1),1-((1-cV)*scJb));
        jbfill(1:Nt,Btint_PVVtau(g,:,2),Btint_PVVtau(g,:,1),1-((1-cVt)*scJb));
        hold on
        h1=plot(1:Nt,K_PVVtau(g,:),'-o','color',cP,'markerfacecolor',cP,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        h2=plot(1:Nt,B_PVVtau(g,:),'-^','color',cV,'markerfacecolor',cV,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        h3=plot(1:Nt,Bt_PVVtau(g,:),'-s','color',cVt,'markerfacecolor',cVt,'markeredgecolor','none','linewidth',lw,'markersize',ms);
        line([1  Nt] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
        ylim([-0.2 1]);
        xlim([1  Nt]);
        text(textloc(1),textloc(2),title,'fontsize' , tfs,'fontweight' , 'b')
        xlabel('Force Channel Trial Number','fontsize' , 20, 'fontweight' , 'b');
        ylabel('Gain','fontsize' , 20, 'fontweight' , 'b');
        l=legend([h1,h2,h3], 'Pos', 'Vel','Delayed Vel');
        set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northeast');
        legend boxoff
        set(gca,'box','off', 'FontSize', afs,'fontweight','b','ytick',0:0.2:1)
        
        figure('Position',fsg)
        hold on
        for i=1:Nt_adapt
            p(i)=plot(K_PVVtau(g,trVec_adapt(i)),B_PVVtau(g,trVec_adapt(i)),'.','markersize',30,'color',c_pert*cm(i),'linewidth',lw);
        end
        text(15,0.9,title,'fontsize' , tfs,'fontweight' , 'b')
        xlabel('Position Gain','fontsize' , 20, 'fontweight' , 'b');
        ylabel('Velocity Gain','fontsize' , 20, 'fontweight' , 'b');
        set(gca,'xtick',0:0.2:0.2,'box','off', 'FontSize', afs,'fontweight','b')
        axis equal
        ylim([-0.1 1]);
        xlim([-0.1 0.3]);
        trColMap=cm*c_pert;
        colormap(trColMap);
        hnd=colorbar;
        ylabel(hnd,'Force Channel Trial Number','fontsize',20);
        set(hnd,'box','off','ticks',[0.2:0.2:1],'yticklabel',5:5:25);
        
        figure('Position',fsg)
        hold on
        for i=1:Nt_adapt
            p(i)=plot(K_PVVtau(g,trVec_adapt(i)),Bt_PVVtau(g,trVec_adapt(i)),'.','markersize',30,'color',c_pert*cm(i),'linewidth',lw);
        end
        text(15,0.9,title,'fontsize' , tfs,'fontweight' , 'b')
        xlabel('Position Gain','fontsize' , 20, 'fontweight' , 'b');
        ylabel('Delayed Velocity Gain','fontsize' , 20, 'fontweight' , 'b');
        set(gca,'xtick',0:0.2:0.2,'box','off', 'FontSize', afs,'fontweight','b')
        axis equal
        ylim([-0.1 1]);
        xlim([-0.1 0.3]);
        trColMap=cm*c_pert;
        colormap(trColMap);
        hnd=colorbar;
        ylabel(hnd,'Force Channel Trial Number','fontsize',20);
        set(hnd,'box','off','ticks',[0.2:0.2:1],'yticklabel',5:5:25);
        
        figure('Position',fs)
        hold on
        for i=1:Nt_adapt
            p(i)=plot(B_PVVtau(g,trVec_adapt(i)),Bt_PVVtau(g,trVec_adapt(i)),'.','markersize',30,'color',c_pert*cm(i),'linewidth',lw);
        end
        text(15,0.9,title,'fontsize' , tfs,'fontweight' , 'b')
        xlabel('Velocity Gain','fontsize' , 20, 'fontweight' , 'b');
        ylabel('Delayed Velocity Gain','fontsize' , 20, 'fontweight' , 'b');
        set(gca,'xtick',0:0.2:1,'box','off', 'FontSize', afs,'fontweight','b')
        axis equal
        ylim([-0.1 1]);
        xlim([-0.1 1]);
        trColMap=cm*c_pert;
        colormap(trColMap);
        hnd=colorbar;
        ylabel(hnd,'Force Channel Trial Number','fontsize',20);
        set(hnd,'box','off','ticks',[0.2:0.2:1],'yticklabel',5:5:25);
        
    end
    
    % Total trials that were not analyzed - both according to duration and
    % multicolinearity
    
    indTrialsRemove=union(find(belowLowCut),find(islargeVIF));
    TrialsRemove=zeros(Ns,Nt);
    TrialsRemove(indTrialsRemove)=1;
    
    percentTrialNotAnalyzed_Group(g)=100*sum(sum(TrialsRemove))/(Nt*Ns);
end

percentTrialNotAnalyzed=100*sum(TrBelowLowCut)/(3*(Nt*Ns)); % percent of trial pairs that were removed from the analysis
              % (out of all FC-BeFC trials in the entire experiment and in all three groups)