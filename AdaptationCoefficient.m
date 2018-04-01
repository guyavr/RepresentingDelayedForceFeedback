function AdaptationCoefficient

clc
clear
close all

c = colormap(lines(10));
close
cnd=c(1,:)*1.2;
cd70=c(3,:)*1.05;
cd100=c(7,:)*1.05;
cd70SF=[255,113,116]/255;

scJb=0.5; % scale color for lightening jbfill

lw=3;
afs=24;
tickfs=20;
lfs=21;

for g=1:4 % 1- ND, 2- D70, 3- D100, 4- D70SF
    
    if g==1
        dat=load('FC_BeFC_AllExpData_ND');
        co=cnd;
    elseif g==2
        dat=load('FC_BeFC_AllExpData_D70');
        co=cd70;
    elseif g==3
        dat=load('FC_BeFC_AllExpData_D100');
        co=cd100;
    else
        dat=load('FC_BeFC_AllExpData_D70SF');
        co=cd70SF;
    end
    
    Ns=size(dat.FC,1); % number of subjects
    
    if g==4
        trVec_adapt=5:29; % force channel trials of only adaptation
    else
        trVec_adapt=6:30; % force channel trials of only adaptation
    end
    Nt_adapt=length(trVec_adapt);
    
    % Find maximum force vector length to create an initialized vector
    TrialLenAdapt_FC=nan(Ns,Nt_adapt);
    TrialLenAdapt_BeFC=nan(Ns,Nt_adapt);
    b_byTrial=NaN(Nt_adapt,2);
    bint_byTrial=NaN(Nt_adapt,2,2);
    b_all=NaN(Ns,Nt_adapt);
    R2_all=NaN(Ns,Nt_adapt);
    
    reg_trialSubj=struct([]);
    
    for i=1:Nt_adapt
        
        for s=1:Ns
            TrialLenAdapt_FC(s,i)=length(dat.FC(s,trVec_adapt(i)).Force);
            TrialLenAdapt_BeFC(s,i)=length(dat.BeFC(s,trVec_adapt(i)).Force);
        end
    end
    
    % Analyse the statistics of the duration of all FC and BeFC trials (for
    % all subjects)
    TrialLenAdapt_All=[reshape(TrialLenAdapt_FC,numel(TrialLenAdapt_FC),1); reshape(TrialLenAdapt_BeFC,numel(TrialLenAdapt_BeFC),1)];
    
    Cut = prctile(TrialLenAdapt_All,10); % for 10% cut off
    lowCut=ceil(Cut(1));
    
    % Create matrices for TrialLen_FC and for TrialLen_BeFC that contain 1
    % for trials shorter than lowCut
    belowLowCut_FC=zeros(Ns,Nt_adapt);
    belowLowCut_BeFC=zeros(Ns,Nt_adapt);
    belowLowCut_FC(TrialLenAdapt_FC<lowCut)=1;
    belowLowCut_BeFC(TrialLenAdapt_BeFC<lowCut)=1;
    
    indLowCut=union(find(belowLowCut_FC),find(belowLowCut_BeFC));
    belowLowCut=zeros(Ns,Nt_adapt);
    belowLowCut(indLowCut)=1;
    
    UniteTrLen=lowCut-1; % It is -1 because of the acceleration vector
    
    TrPairRemove=zeros(Ns,Nt_adapt); % pairs of BeFC-FC trials to remove when the force during BeFC was always zero;
    
    for i=1:Nt_adapt
        % To find a b value of all subjects in one trial- concatenate the
        % forces of all subjects;
        Force_FC_trial=[];
        Force_BeFC_trial=[];
        
        for s=1:Ns
            % calculate the regression with the original length
            TrLen=min(TrialLenAdapt_FC(s,i),TrialLenAdapt_BeFC(s,i));
            
            Force_FC=dat.FC(s,trVec_adapt(i)).Force(1:TrLen)';
            Force_BeFC=dat.BeFC(s,trVec_adapt(i)).Force(1:TrLen)';
            
            FreeVar_all=ones(TrLen,1);
            [reg_trialSubj(s,i).b,~,reg_trialSubj(s,i).r,reg_trialSubj(s,i).rint,reg_trialSubj(s,i).stats]=regress(Force_FC,[Force_BeFC FreeVar_all]);
            b_all(s,i)=reg_trialSubj(s,i).b(1); % Calculate the B value for each trial for each subject (for statistical analysis)
            R2_all(s,i)=reg_trialSubj(s,i).stats(1); % Calculate the B value for each trial for each subject (for statistical analysis)
            
            if isempty(Force_BeFC) % if a perturbing force was not applied
                TrPairRemove(s,i)=1;
            end
            
            if TrLen>=UniteTrLen % keep only movements that are not too fast
                Force_FC_trial=[Force_FC_trial;Force_FC];
                Force_BeFC_trial=[Force_BeFC_trial;Force_BeFC];
            end
            
        end
        
        T=size(Force_FC_trial,1); % Number of elements in the response vector
        FreeVar=ones(T,1);
        [b_byTrial(i,:),bint_byTrial(i,:,:)]=regress(Force_FC_trial,[Force_BeFC_trial FreeVar]);
        
    end
    
    b_all(TrPairRemove==1)=nan;

    NtCon=5; % number of trials to average for each subject
    trEA=1:NtCon; % Early Adaptation trials
    trLA=Nt_adapt-NtCon+1:Nt_adapt; % Late Adaptation trials
    
    b_EA=b_all(:,trEA)'; % transpose to average for each subject
    b_LA=b_all(:,trLA)'; % transpose to average for each subject
    
    mb_EA=nanmean(b_EA)'; % average for each subject
    mb_LA=nanmean(b_LA)'; % average for each subject
    
    if g==1
        b_byTrial_ND=b_byTrial;
        bint_byTrial_ND=bint_byTrial;
        mb_EA_ND=mb_EA;
        mb_LA_ND=mb_LA;
    elseif g==2
        b_byTrial_D70=b_byTrial;
        bint_byTrial_D70=bint_byTrial;
        mb_EA_D70=mb_EA;
        mb_LA_D70=mb_LA;
    elseif g==3
        b_byTrial_D100=b_byTrial;
        bint_byTrial_D100=bint_byTrial;
        mb_EA_D100=mb_EA;
        mb_LA_D100=mb_LA;
    else
        b_byTrial_D70SF=b_byTrial;
        bint_byTrial_D70SF=bint_byTrial;
        mb_EA_D70SF=mb_EA;
        mb_LA_D70SF=mb_LA;
    end
    
end

for g=1:2 % 1- ND, D70 and D100; 2- D70_SF
    if g==1
        Ns=length(mb_EA_ND);
        scrsz = [1,1,1280,800];
        figure('Position',[100 100 scrsz(3)*2/4 scrsz(4)*2.5/4])
        hold on
        p1=plot(1:Nt_adapt, b_byTrial_ND(:,1), 'Color' ,cnd, 'linewidth' ,lw);
        p2=plot(1:Nt_adapt, b_byTrial_D70(:,1), 'Color' ,cd70, 'linewidth' ,lw);
        p3=plot(1:Nt_adapt, b_byTrial_D100(:,1), 'Color' ,cd100, 'linewidth' ,lw);
        jbfill(1:Nt_adapt,bint_byTrial_ND(:,1,2)',bint_byTrial_ND(:,1,1)',1-((1-cnd)*scJb));
        jbfill(1:Nt_adapt,bint_byTrial_D70(:,1,2)',bint_byTrial_D70(:,1,1)',1-((1-cd70)*scJb));
        jbfill(1:Nt_adapt,bint_byTrial_D100(:,1,2)',bint_byTrial_D100(:,1,1)',1-((1-cd100)*scJb));
        ylim([0  .8]);
        xlim([1  Nt_adapt]);
        hold on
        h1=plot(1:Nt_adapt, b_byTrial_ND(:,1), 'Color' ,cnd, 'linewidth' ,lw);
        h2=plot(1:Nt_adapt, b_byTrial_D70(:,1), 'Color' ,cd70, 'linewidth' ,lw);
        h3=plot(1:Nt_adapt, b_byTrial_D100(:,1), 'Color' ,cd100, 'linewidth' ,lw);
        
        xlabel('Force Channel Trial Number','fontsize' , afs, 'fontweight' , 'b');
        ylabel('Adaptation Coefficient','fontsize' , afs, 'fontweight' , 'b');
        set(gca,'xtick',[5:5:25],'ytick',[0:0.2:0.8],'box','off', 'FontSize', tickfs,'fontweight','b')
        l=legend ([p1,p2,p3],'Group ND','Group D70','Group D100','location','southeast');
        set(l, 'Box', 'off','fontsize',lfs,'fontweight' , 'b')
        
        mmb_EA_ND=nanmean(mb_EA_ND);
        mmb_LA_ND=nanmean(mb_LA_ND);
        mmb_EA_D70=nanmean(mb_EA_D70);
        mmb_LA_D70=nanmean(mb_LA_D70);
        mmb_EA_D100=nanmean(mb_EA_D100);
        mmb_LA_D100=nanmean(mb_LA_D100);
        mm_ND=[mmb_EA_ND,mmb_LA_ND];
        mm_D70=[mmb_EA_D70,mmb_LA_D70];
        mm_D100=[mmb_EA_D100,mmb_LA_D100];
        
        mm=[mm_ND;mm_D70;mm_D100];
        
        % confidence interval (the distance from the mean)
        cimb_EA_ND=nanstd(mb_EA_ND)*1.96/sqrt(Ns);
        cimb_LA_ND=nanstd(mb_LA_ND)*1.96/sqrt(Ns);
        cimb_EA_D70=nanstd(mb_EA_D70)*1.96/sqrt(Ns);
        cimb_LA_D70=nanstd(mb_LA_D70)*1.96/sqrt(Ns);
        cimb_EA_D100=nanstd(mb_EA_D100)*1.96/sqrt(Ns);
        cimb_LA_D100=nanstd(mb_LA_D100)*1.96/sqrt(Ns);
        cim_ND=[cimb_EA_ND,cimb_LA_ND];
        cim_D70=[cimb_EA_D70,cimb_LA_D70];
        cim_D100=[cimb_EA_D100,cimb_LA_D100];
        
        cim=[cim_ND;cim_D70;cim_D100];
        
        plotStatForceRegAdapt_ND_D70_D100(mm,cim)
    else
        scrsz = [1,1,1280,800];
        figure('Position',[100 100 scrsz(3)*2/4 scrsz(4)*2.5/4])
        hold on
        p1=plot(1:Nt_adapt, b_byTrial_D70SF(:,1), 'Color' ,cd70SF, 'linewidth' ,lw);
        jbfill(1:Nt_adapt,bint_byTrial_D70SF(:,1,2)',bint_byTrial_D70SF(:,1,1)',1-((1-cd70SF)*scJb));
        ylim([0  .8]);
        % ylim([0  1]);
        xlim([1  Nt_adapt]);
        hold on
        h1=plot(1:Nt_adapt, b_byTrial_D70SF(:,1), 'Color' ,cd70SF, 'linewidth' ,lw);
        
        xlabel('Force Channel Trial Number','fontsize' , afs-4, 'fontweight' , 'b');
        ylabel('Adaptation Coefficient','fontsize' , afs-4, 'fontweight' , 'b');
        set(gca,'xtick',[5:5:25],'ytick',[0:0.2:0.8],'box','off', 'FontSize', tickfs,'fontweight','b')
        l=legend (p1,'Group D70\_SF','location','southeast');
        set(l, 'Box', 'off','fontsize',lfs,'fontweight' , 'b')
    end
end

