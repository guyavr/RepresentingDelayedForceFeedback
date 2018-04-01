function RepresentationAnalysis_LateAdaptation

clc
clear
close all

c = colormap(lines(10));
close
cnd_pert=c(1,:)*1.2;
cd70_pert=c(3,:)*1.05;
cd100_pert=c(7,:)*1.05;
cd70SF_pert=[255,113,116]/255;

dF_rep=0.4; % darkening factor
dF_act=0.65; % darkening factor
scJb=0.5; % scale color for lightening jbfill

cP=c(2,:)*1.167;
cV=c(5,:)*0.77;
cA=c(4,:)*0.9;
cVt=c(6,:)*0.4;

% colors of individuals' actual forces
cnd_subj=[0 178 238;0 104 139; 28 134 238; 16 78 139; 0 0 205; 65 0105 225; 39 64 139; 79 148 205; 51 102 153; 3 180 204]/255;
cd70_subj=[184 134 11; 255 215 0; 238 201 0; 230 210 120; 255 127 0; 205 133 0; 160 82 45; 255 200 40; 255 185 15; 139 90 0]/255;
cd100_subj=[205 16 118; 160 0 0; 205 96 50; 205 85 85; 225 48 48; 255 99 71; 255 106 106; 255 62 150; 255 0 0; 255 121 156]/255;
cd70SF_subj=[255,20,147; 250,128,114; 205,92,92; 219,112,147; 220,20,60; 240,128,128; 233,150,122; 255,113,116]/255;
s2ms=1000;

afs=19;
tfs=20;
lfs=18;

ofs = get(gcf,'Position'); % original figure size
close
sfs=1.3; % scale figure size
mfs=3; % move figure
fs=[ofs(1:2)/mfs ofs(3:4)*sfs];% figure size
fs_bar3=[ofs(1:2)/mfs ofs(3)*sfs/2 ofs(4)*sfs];% figure size
fs_hist_nPeaks=[ofs(1:2)/mfs ofs(3)*sfs/2.5 ofs(4)*sfs/1.2];% figure size
fs_hist_nPeaksTimes=[ofs(1:2)/mfs ofs(3)*sfs/1.2 ofs(4)*sfs/1.2];% figure size
    
lwRep=3;
lw=4;

xLim=[0 0.545]*s2ms;
textloc=[0.35*s2ms 3.3]; % location of text box

ytick=[0:1:3];
xtick=[0:100:500];
xticklabel=[-100:100:400];

leg_act='Actual Force';
leg_pert='Perturbation';

st=0.005;

[bwb,bwa]=butter(2,0.1); % build lowpass butterworth filter for the acceleration signal- 2nd order for 10Hz (10Hz/200Hz*2)

trAna=16:25; % end of adaptation trials - for representation analysis.
Nt_ana=length(trAna);

% Changing the data into force units by multiplying with the
% appropriate factors
VFactor=0.06; % N/cm

maxForceBase_FC=cell(1,4);
mMaxForceBase_FC=cell(1,4);
mmMaxForceBase_FC=nan(1,4);

nPks_ForceFC=cell(1,4);
timepks=cell(1,4);
timepksVec=cell(1,4);
dTimepks=cell(1,4);
dTimepksVec=cell(1,4);

maxVel=cell(1,4);
maxAcc=cell(1,4);
mMaxVel_lateAdapt=cell(1,4);

for g=1:4 % 1- ND, 2- D70, 3- D100, 4- D70SF
    
    if g==1
        dat=load('FC_BeFC_AllExpData_ND');
        c_pert=cnd_pert;
        c_subj=cnd_subj;
        delVal=0;
        title='Group ND';
        PFactor=0.38; % N*s/cm
        AFactor=0.0056; % N*s^2/cm
    elseif g==2
        dat=load('FC_BeFC_AllExpData_D70');
        c_pert=cd70_pert;
        c_subj=cd70_subj;
        delVal=0.07;
        title='Group D70';
        PFactor=0.32; % N*s/cm
        AFactor=0.0068; % N*s^2/cm
    elseif g==3
        dat=load('FC_BeFC_AllExpData_D100');
        c_pert=cd100_pert;
        c_subj=cd100_subj;
        delVal=0.1;
        title='Group D100';
        PFactor=0.26; % N*s/cm
        AFactor=0.0073; % N*s^2/cm
    else
        dat=load('FC_BeFC_AllExpData_D70SF');
        c_pert=cd70SF_pert;
        c_subj=cd70SF_subj;
        delVal=0.07;
        title='Group D70\_SF';
        PFactor=0.21; % N*s/cm
        AFactor=0.0082; % N*s^2/cm
        textloc=[300 2.6];
    end
    
    if g==4
        trVec_base=1:4;
        trVec_adapt=5:29; % force channel trials of only adaptation
        yLim=[-0.6 2.8];
    else
        trVec_base=1:5;
        trVec_adapt=6:30; % force channel trials of only adaptation
        yLim=[-0.5 3.6];
    end
    
    Nt_base=length(trVec_base);
    Nt_adapt=length(trVec_adapt);

    delSamp=floor(delVal/st);
    
    c_act=c_pert*dF_act;
    c_rep=c_pert*dF_rep;
    
    Ns=size(dat.FC,1); % number of subjects
    
    % Arrange actual force data according to trial- for each trial an array of
    % Ns (number of subjects) force profiles (in raws) - only for adaptation
    
    % Find maximum force vector length to create an initialized vector
    TrialLenAdapt_FC=nan(Ns,Nt_adapt);
    TrialLenAdapt_BeFC=nan(Ns,Nt_adapt);
    
    % Calculate max force at force channel during baseline (to find the
    % minimum applied force when no representation is there - for cutoff in
    % findpeaks)
    maxForceBase_FC{g}=nan(Ns,Nt_base);
    for i=1:Nt_base
        for s=1:Ns
            maxForceBase_FC{g}(s,i)=max(dat.FC(s,i).Force);
        end
    end
    mMaxForceBase_FC{g}=mean(maxForceBase_FC{g},2);
    mmMaxForceBase_FC(g)=mean(mMaxForceBase_FC{g});
    
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
    
    UniteTrLen(g)=lowCut-1; % It is -1 because of the acceleration vector
    
    TrPairRemove=zeros(Ns,Nt_adapt); % pairs of BeFC-FC trials to remove when the force during BeFC was always zero;
    
    maxVel{g}=nan(Ns,Nt_adapt); % find max velocity for VFactor
    maxAcc{g}=nan(Ns,Nt_adapt); % find max velocity for AFactor
    
    for i=1:Nt_adapt
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
            % calculate the regression with the original length 
            TrLen=min(TrialLenAdapt_FC(s,i),TrialLenAdapt_BeFC(s,i))-1;
            
            Force_FC=dat.FC(s,trVec_adapt(i)).Force(1:TrLen)';
            
            Force_BeFC=dat.BeFC(s,trVec_adapt(i)).Force(1:TrLen)';
            PosY=dat.BeFC(s,trVec_adapt(i)).PosY(1:TrLen)';
            VelY=dat.BeFC(s,trVec_adapt(i)).VelY(1:TrLen)';
            AccY_org=diff(dat.BeFC(s,trVec_adapt(i)).VelY(1:TrLen+1)')/st; % acceleration before filtering
            AccY=filtfilt(bwb,bwa,AccY_org); % acceleration after filtering
            VeltauY=[zeros(1,delSamp) dat.BeFC(s,trVec_adapt(i)).VelY(1:TrLen-delSamp)]';
            if ~Force_BeFC % if a perturbing force was not applied
                TrPairRemove(s,i)=1;
            end
            
            if Force_BeFC
                % Calculate for only BeFC trials
                maxVel{g}(s,i)=max(VelY);
                maxAcc{g}(s,i)=max(abs(AccY));
            end
            
            % for calculating the primitives's weight for each trial
            if TrLen>=UniteTrLen(g) % keep only movements that are not too fast
                Force_FC_trial=[Force_FC_trial;Force_FC];
                Force_BeFC_trial=[Force_BeFC_trial;Force_BeFC];
                PosY_trial=[PosY_trial;PosY];
                VelY_trial=[VelY_trial;VelY];
                AccY_trial=[AccY_trial;AccY];
                VelYtau_trial=[VelYtau_trial;VeltauY];
                
                % To create the primitives matrix- should be with dummy variables
                % according to the subjects analyzed in each trial
                trLenAnalyzeTrial=[trLenAnalyzeTrial,TrLen]; % subjects to analyze is length(trLenAnalyzeTrial)
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
        
        % vel (V)
        prim=[V]; primdum=[Vdum];
        [R2_V(g,i),LogL_V,Npar_V,B_V(g,i),Bint_V(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % delayed vel (Vtau)
        prim=[Vtau]; primdum=[Vtaudum];
        [R2_Vtau(g,i),LogL_Vtau,Npar_Vtau,Bt_Vtau(g,i),Btint_Vtau(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % pos (P) and vel (V)
        prim=[P V]; primdum=[Pdum Vdum];
        [R2_PV(g,i),LogL_PV,Npar_PV,K_PV(g,i),Kint_PV(g,i,:),B_PV(g,i),Bint_PV(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % pos (P) and acc (A)
        prim=[P A]; primdum=[Pdum Adum];
        [R2_PA(g,i),LogL_PA,Npar_PA,K_PA(g,i),Kint_PA(g,i,:),M_PA(g,i),Mint_PA(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % pos (P) and delayed vel (Vtau)
        prim=[P Vtau]; primdum=[Pdum Vtaudum];
        [R2_PVtau(g,i),LogL_PVtau,Npar_PVtau,K_PVtau(g,i),Kint_PVtau(g,i,:),Bt_PVtau(g,i),Btint_PVtau(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % vel (V) and acc (A)
        prim=[V A]; primdum=[Vdum Adum];
        [R2_VA(g,i),LogL_VA,Npar_VA,B_VA(g,i),Bint_VA(g,i,:),M_VA(g,i),Mint_VA(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % vel (V) and delayed vel (Vtau)
        prim=[V Vtau]; primdum=[Vdum Vtaudum];
        [R2_VVtau(g,i),LogL_VVtau,Npar_VVtau,B_VVtau(g,i),Bint_VVtau(g,i,:),Bt_VVtau(g,i),Btint_VVtau(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % pos (P), vel (V) and acc (A)
        prim=[P V A]; primdum=[Pdum Vdum Adum];
        [R2_PVA(g,i),LogL_PVA,Npar_PVA,K_PVA(g,i),Kint_PVA(g,i,:),B_PVA(g,i),Bint_PVA(g,i,:),M_PVA(g,i),Mint_PVA(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % pos (P), vel (V) and delayed vel (Vtau)
        prim=[P V Vtau]; primdum=[Pdum Vdum Vtaudum];
        [R2_PVVtau(g,i),LogL_PVVtau,Npar_PVVtau,K_PVVtau(g,i),Kint_PVVtau(g,i,:),B_PVVtau(g,i),Bint_PVVtau(g,i,:),Bt_PVVtau(g,i),Btint_PVVtau(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % pos (P), vel (V), delayed vel (Vtau) and acc (A)
        prim=[P V Vtau A]; primdum=[Pdum Vdum Vtaudum Adum];
        [R2_PVVtauA(g,i),LogL_PVVtauA,Npar_PVVtauA,K_PVVtauA(g,i),Kint_PVVtauA(g,i,:),B_PVVtauA(g,i),Bint_PVVtauA(g,i,:),Bt_PVVtauA(g,i),Btint_PVVtauA(g,i,:),M_PVVtauA(g,i),Mint_PVVtauA(g,i,:)] = forcePrimReg(FAct,prim,primdum);
        
        % LogL for all models
        LogL=[LogL_V LogL_Vtau LogL_PV LogL_PA LogL_PVtau LogL_VA LogL_VVtau LogL_PVA LogL_PVVtau LogL_PVVtauA];
        Npar=[Npar_V Npar_Vtau Npar_PV Npar_PA Npar_PVtau Npar_VA Npar_VVtau Npar_PVA Npar_PVVtau Npar_PVVtauA];
        
        T=length(FAct)*ones(1,length(Npar));
        
        [aic(g,i,:),bic(g,i,:)]=aicbic(LogL,Npar,T);
        
    end
    
    mMaxVel_lateAdapt{g}=nanmean(maxVel{g},2);
    
    % to calculate the entire weight for the end of adaptation, I
    % first need to concatenate for each subject each primitive
    % from the trAna trials, and then concatenate the subjects
    Force_FC_lateAdapt=[];
    Force_BeFC_lateAdapt=[];
    PosY_lateAdapt=[];
    VelY_lateAdapt=[];
    AccY_lateAdapt=[];
    VelYtau_lateAdapt=[];
    
    trLenAnalyzeTrial_lateAdapt=[];
    
    nPks_ForceFC{g}=nan(Ns,Nt_ana);
    timepks{g}=cell(Ns,Nt_ana);
    timepksVec{g}=[];
    dTimepks{g}=cell(Ns,Nt_ana);
    dTimepksVec{g}=[];
    
    % Find the properties- force and velocity max time and offset
    MaxActForce=nan(Ns,Nt_ana);
    TMaxActForce=nan(Ns,Nt_ana);
    TStartActForce=nan(Ns,Nt_ana);
    MaxVel=nan(Ns,Nt_ana);
    TMaxVel=nan(Ns,Nt_ana);
    TStartVel=nan(Ns,Nt_ana);
    MaxVeltau=nan(Ns,Nt_ana);
    TMaxVeltau=nan(Ns,Nt_ana);
    TStartVeltau=nan(Ns,Nt_ana);
    
    for s=1:Ns
        
        Force_FC_lateAdapt_subj=[];
        Force_BeFC_lateAdapt_subj=[];
        PosY_lateAdapt_subj=[];
        VelY_lateAdapt_subj=[];
        AccY_lateAdapt_subj=[];
        VelYtau_lateAdapt_subj=[];
        
        for i=1:Nt_ana
            
            % calculate the regression with the original length 
            TrLen=min(TrialLenAdapt_FC(s,trAna(i)),TrialLenAdapt_BeFC(s,trAna(i)))-1;
            
            t=st*(1:TrLen)*s2ms;
    
            Force_FC=dat.FC(s,trVec_adapt(trAna(i))).Force(1:TrLen)';
            
            Force_BeFC=dat.BeFC(s,trVec_adapt(trAna(i))).Force(1:TrLen)';
            PosY=dat.BeFC(s,trVec_adapt(trAna(i))).PosY(1:TrLen)';
            VelY=dat.BeFC(s,trVec_adapt(trAna(i))).VelY(1:TrLen)';
            AccY_org=diff(dat.BeFC(s,trVec_adapt(trAna(i))).VelY(1:TrLen+1)')/st; % acceleration before filtering
            AccY=filtfilt(bwb,bwa,AccY_org); % acceleration after filtering
            VeltauY=[zeros(1,delSamp) dat.BeFC(s,trVec_adapt(trAna(i))).VelY(1:TrLen-delSamp)]';
            if ~Force_BeFC % if a perturbing force was not applied
                TrPairRemove(s,i)=1;
            end
            
            % for calculating the primitives's weight for each trial
            if TrLen>=UniteTrLen(g) % keep only movements that are not too fast
                filtForce_FC=filtfilt(bwb,bwa,Force_FC);
                [pks,pksloc] = findpeaks(filtForce_FC,'MinPeakHeight',mMaxForceBase_FC{g}(s));
                nPks_ForceFC{g}(s,i)=length(pks);
                timepks{g}{s,i}=pksloc*st*s2ms;
                timepksVec{g}=[timepksVec{g} timepks{g}{s,i}'];
                dTimepks{g}{s,i}=diff(pksloc)*st*s2ms;
                dTimepksVec{g}=[dTimepksVec{g} dTimepks{g}{s,i}'];
                
                % Find the properties- force and velocity max time and offset
                [MaxActForce(s,i),TMaxActForce(s,i)]=findMaxProfile(t,Force_FC);
                TStartActForce(s,i)=findStartProfile(t,Force_FC);
                [MaxVel(s,i),TMaxVel(s,i)]=findMaxProfile(t,VelY);
                TStartVel(s,i)=findStartProfile(t,VelY);
                [MaxVeltau(s,i),TMaxVeltau(s,i)]=findMaxProfile(t,VeltauY);
                TStartVeltau(s,i)=findStartProfile(t,VeltauY);
                
                
                Force_FC_lateAdapt_subj=[Force_FC_lateAdapt_subj;Force_FC];
                Force_BeFC_lateAdapt_subj=[Force_BeFC_lateAdapt_subj;Force_BeFC];
                PosY_lateAdapt_subj=[PosY_lateAdapt_subj;PosY];
                VelY_lateAdapt_subj=[VelY_lateAdapt_subj;VelY];
                AccY_lateAdapt_subj=[AccY_lateAdapt_subj;AccY];
                VelYtau_lateAdapt_subj=[VelYtau_lateAdapt_subj;VeltauY];
                
            end
            
        end
        
        Force_FC_lateAdapt=[Force_FC_lateAdapt;Force_FC_lateAdapt_subj];
        Force_BeFC_lateAdapt=[Force_BeFC_lateAdapt;Force_BeFC_lateAdapt_subj];
        PosY_lateAdapt=[PosY_lateAdapt;PosY_lateAdapt_subj];
        VelY_lateAdapt=[VelY_lateAdapt;VelY_lateAdapt_subj];
        AccY_lateAdapt=[AccY_lateAdapt;AccY_lateAdapt_subj];
        VelYtau_lateAdapt=[VelYtau_lateAdapt;VelYtau_lateAdapt_subj];
        
        lateAdaptLen=length(Force_FC_lateAdapt_subj);
        % To create the primitives matrix with dummy variables
        % according to the subjects analyzed in each trial
        trLenAnalyzeTrial_lateAdapt=[trLenAnalyzeTrial_lateAdapt,lateAdaptLen]; 
        
    end
    
    % Average ten last trials for each subject
    mMaxActForce=nanmean(MaxActForce');
    mTMaxActForce=nanmean(TMaxActForce');
    mTStartActForce=nanmean(TStartActForce');
    mMaxVel=nanmean(MaxVel');
    mTMaxVel=nanmean(TMaxVel');
    mTStartVel=nanmean(TStartVel');
    mMaxVeltau=nanmean(MaxVeltau');
    mTMaxVeltau=nanmean(TMaxVeltau');
    mTStartVeltau=nanmean(TStartVeltau');
    
    FAct=Force_FC_lateAdapt;
    P=PosY_lateAdapt*PFactor;
    V=VelY_lateAdapt*VFactor;
    A=AccY_lateAdapt*AFactor;
    Vtau=VelYtau_lateAdapt*VFactor;
    
    Pdum = primDummy_differentTrLen(P,trLenAnalyzeTrial_lateAdapt);
    Vdum = primDummy_differentTrLen(V,trLenAnalyzeTrial_lateAdapt);
    Vtaudum = primDummy_differentTrLen(Vtau,trLenAnalyzeTrial_lateAdapt);
    Adum = primDummy_differentTrLen(A,trLenAnalyzeTrial_lateAdapt);
    
    % vel (V)
    prim=[V]; primdum=[Vdum];
    [R2_V_lateAdapt(g),LogL_V_lateAdapt,Npar_V_lateAdapt,B_V_lateAdapt(g),Bint_V_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Bs_V_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % delayed vel (Vtau)
    prim=[Vtau]; primdum=[Vtaudum];
    [R2_Vtau_lateAdapt(g),LogL_Vtau_lateAdapt,Npar_Vtau_lateAdapt,Bt_Vtau_lateAdapt(g),Btint_Vtau_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Bts_Vtau_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % pos (P) and vel (V)
    prim=[P V]; primdum=[Pdum Vdum];
    [R2_PV_lateAdapt(g),LogL_PV_lateAdapt,Npar_PV_lateAdapt,K_PV_lateAdapt(g),Kint_PV_lateAdapt(g,:),B_PV_lateAdapt(g),Bint_PV_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Ks_PV_lateAdapt{g},Bs_PV_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % pos (P) and acc (A)
    prim=[P A]; primdum=[Pdum Adum];
    [R2_PA_lateAdapt(g),LogL_PA_lateAdapt,Npar_PA_lateAdapt,K_PA_lateAdapt(g),Kint_PA_lateAdapt(g,:),M_PA_lateAdapt(g),Mint_PA_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Ks_PA_lateAdapt{g},Ms_PA_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % pos (P) and delayed vel (Vtau)
    prim=[P Vtau]; primdum=[Pdum Vtaudum];
    [R2_PVtau_lateAdapt(g),LogL_PVtau_lateAdapt,Npar_PVtau_lateAdapt,K_PVtau_lateAdapt(g),Kint_PVtau_lateAdapt(g,:),Bt_PVtau_lateAdapt(g),Btint_PVtau_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Ks_PVtau_lateAdapt{g},Bts_PVtau_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % vel (V) and acc (A)
    prim=[V A]; primdum=[Vdum Adum];
    [R2_VA_lateAdapt(g),LogL_VA_lateAdapt,Npar_VA_lateAdapt,B_VA_lateAdapt(g),Bint_VA_lateAdapt(g,:),M_VA_lateAdapt(g),Mint_VA_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Bs_VA_lateAdapt{g},Ms_VA_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % vel (V) and delayed vel (Vtau)
    prim=[V Vtau]; primdum=[Vdum Vtaudum];
    [R2_VVtau_lateAdapt(g),LogL_VVtau_lateAdapt,Npar_VVtau_lateAdapt,B_VVtau_lateAdapt(g),Bint_VVtau_lateAdapt(g,:),Bt_VVtau_lateAdapt(g),Btint_VVtau_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Bs_VVtau_lateAdapt{g},Bts_VVtau_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % pos (P), vel (V) and acc (A)
    prim=[P V A]; primdum=[Pdum Vdum Adum];
    [R2_PVA_lateAdapt(g),LogL_PVA_lateAdapt,Npar_PVA_lateAdapt,K_PVA_lateAdapt(g),Kint_PVA_lateAdapt(g,:),B_PVA_lateAdapt(g),Bint_PVA_lateAdapt(g,:),M_PVA_lateAdapt(g),Mint_PVA_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Ks_PVA_lateAdapt{g},Bs_PVA_lateAdapt{g},Ms_PVA_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % pos (P), vel (V) and delayed vel (Vtau)
    prim=[P V Vtau]; primdum=[Pdum Vdum Vtaudum];
    [R2_PVVtau_lateAdapt(g),LogL_PVVtau_lateAdapt,Npar_PVVtau_lateAdapt,K_PVVtau_lateAdapt(g),Kint_PVVtau_lateAdapt(g,:),B_PVVtau_lateAdapt(g),Bint_PVVtau_lateAdapt(g,:),Bt_PVVtau_lateAdapt(g),Btint_PVVtau_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Ks_PVVtau_lateAdapt{g},Bs_PVVtau_lateAdapt{g},Bts_PVVtau_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % pos (P), vel (V), delayed vel (Vtau) and acc (A)
    prim=[P V Vtau A]; primdum=[Pdum Vdum Vtaudum Adum];
    [R2_PVVtauA_lateAdapt(g),LogL_PVVtauA_lateAdapt,Npar_PVVtauA_lateAdapt,K_PVVtauA_lateAdapt(g),Kint_PVVtauA_lateAdapt(g,:),B_PVVtauA_lateAdapt(g),Bint_PVVtauA_lateAdapt(g,:),Bt_PVVtauA_lateAdapt(g),Btint_PVVtauA_lateAdapt(g,:),M_PVVtauA_lateAdapt(g),Mint_PVVtauA_lateAdapt(g,:)] = forcePrimReg(FAct,prim,primdum);
    [Ks_PVVtauA_lateAdapt{g},Bs_PVVtauA_lateAdapt{g},Bts_PVVtauA_lateAdapt{g},Ms_PVVtauA_lateAdapt{g}] = forcePrimRegIndivPartic(FAct,prim,primdum,Ns);
    
    % LogL for all models
    LogL_lateAdapt=[LogL_V_lateAdapt LogL_Vtau_lateAdapt LogL_PV_lateAdapt LogL_PA_lateAdapt LogL_PVtau_lateAdapt LogL_VA_lateAdapt LogL_VVtau_lateAdapt LogL_PVA_lateAdapt LogL_PVVtau_lateAdapt LogL_PVVtauA_lateAdapt];
    Npar_lateAdapt=[Npar_V_lateAdapt Npar_Vtau_lateAdapt Npar_PV_lateAdapt Npar_PA_lateAdapt Npar_PVtau_lateAdapt Npar_VA_lateAdapt Npar_VVtau_lateAdapt Npar_PVA_lateAdapt Npar_PVVtau_lateAdapt Npar_PVVtauA_lateAdapt];
    
    T_lateAdapt=length(FAct)*ones(1,length(Npar_lateAdapt));
    
    [aic_lateAdapt(g,:),bic_lateAdapt(g,:)]=aicbic(LogL_lateAdapt,Npar_lateAdapt,T_lateAdapt);
    
    % Figures
    t=st*(1:UniteTrLen(g))*s2ms;
    
    % Calculate the mean actual force of the last "trAna" Force Channel
    % trials and the mean perturbation force of the last "trAna" trials
    % preceding the ForceChannel trials. Also, arrange the primitives to
    % plot.
    ActF.all=NaN(Ns*Nt_ana,UniteTrLen(g)); % take all the trAna trials of all subjects together
    PertF.all=NaN(Ns*Nt_ana,UniteTrLen(g)); % take all the trAna trials of all subjects together
    PosY.all=NaN(Ns*Nt_ana,UniteTrLen(g)); % take all the trAna trials of all subjects together
    VelY.all=NaN(Ns*Nt_ana,UniteTrLen(g)); % take all the trAna trials of all subjects together
    VeltauY.all=NaN(Ns*Nt_ana,UniteTrLen(g)); % take all the trAna trials of all subjects together
    AccY.all=NaN(Ns*Nt_ana,UniteTrLen(g)); % take all the trAna trials of all subjects together
    
    ActF.mSubj=NaN(Ns,UniteTrLen(g)); % mean late adaptation actual force for each subject
    
    for s=1:Ns
        for i=1:length(trAna)
            if ~belowLowCut(s,i+trAna(1)-1)
                ActF.all((s-1)*Nt_ana+i,:)=dat.FC(s,trVec_adapt(trAna(i))).Force(1:UniteTrLen(g));
                PertF.all((s-1)*Nt_ana+i,:)=dat.BeFC(s,trVec_adapt(trAna(i))).Force(1:UniteTrLen(g));
                PosY.all((s-1)*Nt_ana+i,:)=dat.BeFC(s,trVec_adapt(trAna(i))).PosY(1:UniteTrLen(g));
                VelY.all((s-1)*Nt_ana+i,:)=dat.BeFC(s,trVec_adapt(trAna(i))).VelY(1:UniteTrLen(g));
                AccY_org=diff(dat.BeFC(s,trVec_adapt(trAna(i))).VelY(1:UniteTrLen(g)+1))/st; % acceleration before filtering
                AccY.all((s-1)*Nt_ana+i,:)=filtfilt(bwb,bwa,AccY_org); % acceleration after filtering
                VeltauY.all((s-1)*Nt_ana+i,:)=[zeros(1,delSamp) dat.BeFC(s,trVec_adapt(trAna(i))).VelY(1:UniteTrLen(g)-delSamp)];
            end
        end
        
        ActF.mSubj(s,:)=nanmean(ActF.all((s-1)*Nt_ana+1:s*Nt_ana,:)); 
    end
    
    nHist{g}=histcounts(nPks_ForceFC{g});
    pHist{g}=nHist{g}/(Ns*Nt_ana); % probability
    
    figure('Position',fs_hist_nPeaks)
    hi=bar(1:length(pHist{g}),pHist{g}, 'facecolor' ,c_pert, 'edgecolor' ,c_pert);
    ylim([0 .55]);
    xlim([0.5 5.5]);
    xlabel('# Peaks','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Probability','fontsize' , 20, 'fontweight' , 'b');
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',1:5,'ytick',0:.1:.5,'layer','top')
    
    figure('Position',fs_hist_nPeaksTimes)
    nHistTime{g}=histcounts(timepksVec{g},0:25:700);
    pHistTime{g}=nHistTime{g}/sum(nHistTime{g}); % probability
    hiTime=histogram(timepksVec{g},0:25:700,'normalization','probability', 'facecolor' ,c_act, 'edgecolor' ,[1 1 1],'facealpha',1);
    xlim([0 550]);
    ylim([0 .2]);
    xlabel('Times of Peak Forces [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Probability','fontsize' , 20, 'fontweight' , 'b');
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',0:.04:.3,'layer','top')
    
    trForSTE=sum(~isnan(ActF.all(:,1))); % number of trials used (without NaN for calculating STE)
    
    ActF.m=nanmean(ActF.all); % mean of all subjects
    ActF.ste=nanstd(ActF.all)/sqrt(trForSTE);
    ActF.ci=NonNanDataCI(ActF.all);
    
    PertF.m=nanmean(PertF.all); % mean of all subjects
    PertF.ste=nanstd(PertF.all)/sqrt(trForSTE);
    PertF.ci=NonNanDataCI(PertF.all);
    
    % Plot the ActForce together with the PerturbationForce
    figure('Position',fs)
    hold on
    p1=plot(t, PertF.m, 'color' ,c_pert  , 'linewidth' ,lw);
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    jbfill(t,PertF.ci(2,:),PertF.ci(1,:),1-((1-c_pert)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    hold on
    h1=plot(t, PertF.m, 'color' ,c_pert  , 'linewidth' ,lw);
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h1,h2], leg_pert ,leg_act,  'Pos','Vel');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    text(textloc(1),textloc(2),title,'fontsize' , tfs,'fontweight' , 'b')
    
    % Plot the ActForce for each subject
    figure('Position',fs)
    hold on
    for s=1:Ns
        p(s)=plot(t, ActF.mSubj(s,:),'-','color',c_subj(s,:), 'linewidth' ,lw); % data
    end
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    text(20,textloc(2),'Actual Force','fontsize' , lfs,'fontweight' , 'b')
    
    % PV
    PredF.PV.KP.all=PFactor*K_PV_lateAdapt(g)*PosY.all;
    PredF.PV.BV.all=VFactor*B_PV_lateAdapt(g)*VelY.all;
    PredF.PV.Tot.all=PredF.PV.KP.all + PredF.PV.BV.all;
    PredF.PV.KP.m=nanmean(PredF.PV.KP.all); % mean of all subjects
    PredF.PV.KP.ste=nanstd(PredF.PV.KP.all)/sqrt(trForSTE);
    PredF.PV.KP.ci=NonNanDataCI(PredF.PV.KP.all);
    PredF.PV.BV.m=nanmean(PredF.PV.BV.all); % mean of all subjects
    PredF.PV.BV.ste=nanstd(PredF.PV.BV.all)/sqrt(trForSTE);
    PredF.PV.BV.ci=NonNanDataCI(PredF.PV.BV.all);
    PredF.PV.Tot.m=nanmean(PredF.PV.Tot.all); % mean of all subjects
    PredF.PV.Tot.ste=nanstd(PredF.PV.Tot.all)/sqrt(trForSTE);
    PredF.PV.Tot.ci=NonNanDataCI(PredF.PV.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p3=plot(t, PredF.PV.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    p4=plot(t, PredF.PV.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    p1=plot(t, PredF.PV.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,PredF.PV.Tot.ci(2,:),PredF.PV.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,PredF.PV.KP.ci(2,:),PredF.PV.KP.ci(1,:),1-((1-cP)*scJb));
    jbfill(t,PredF.PV.BV.ci(2,:),PredF.PV.BV.ci(1,:),1-((1-cV)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h3=plot(t, PredF.PV.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    h4=plot(t, PredF.PV.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    h1=plot(t, PredF.PV.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h3,h4], leg_act, 'Representation Model' , 'Pos','Vel');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    
    % arrange for plotErrorBarGain
    g_PV_lateAdapt=[K_PV_lateAdapt(g) B_PV_lateAdapt(g)];
    gint_PV_lateAdapt=[Kint_PV_lateAdapt(g,:)' Bint_PV_lateAdapt(g,:)'];
    gs_PV_lateAdapt=[Ks_PV_lateAdapt{g} Bs_PV_lateAdapt{g}];
    c_PV=[cP;cV];
    legGain={'Pos','Vel'};
    
    figure('Position',fs_bar3)
    hold on
    plotErrorBarGain(g_PV_lateAdapt,gint_PV_lateAdapt,gs_PV_lateAdapt,c_PV,legGain)
    
    % PVtau
    PredF.PVtau.KP.all=PFactor*K_PVtau_lateAdapt(g)*PosY.all;
    PredF.PVtau.BtVtau.all=VFactor*Bt_PVtau_lateAdapt(g)*VeltauY.all;
    PredF.PVtau.Tot.all=PredF.PVtau.KP.all + PredF.PVtau.BtVtau.all;
    PredF.PVtau.KP.m=nanmean(PredF.PVtau.KP.all); % mean of all subjects
    PredF.PVtau.KP.ste=nanstd(PredF.PVtau.KP.all)/sqrt(trForSTE);
    PredF.PVtau.KP.ci=NonNanDataCI(PredF.PVtau.KP.all);
    PredF.PVtau.BtVtau.m=nanmean(PredF.PVtau.BtVtau.all); % mean of all subjects
    PredF.PVtau.BtVtau.ste=nanstd(PredF.PVtau.BtVtau.all)/sqrt(trForSTE);
    PredF.PVtau.BtVtau.ci=NonNanDataCI(PredF.PVtau.BtVtau.all);
    PredF.PVtau.Tot.m=nanmean(PredF.PVtau.Tot.all); % mean of all subjects
    PredF.PVtau.Tot.ste=nanstd(PredF.PVtau.Tot.all)/sqrt(trForSTE);
    PredF.PVtau.Tot.ci=NonNanDataCI(PredF.PVtau.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p3=plot(t, PredF.PVtau.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    p4=plot(t, PredF.PVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    p1=plot(t, PredF.PVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,PredF.PVtau.Tot.ci(2,:),PredF.PVtau.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,PredF.PVtau.KP.ci(2,:),PredF.PVtau.KP.ci(1,:),1-((1-cP)*scJb));
    jbfill(t,PredF.PVtau.BtVtau.ci(2,:),PredF.PVtau.BtVtau.ci(1,:),1-((1-cVt)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h3=plot(t, PredF.PVtau.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    h4=plot(t, PredF.PVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    h1=plot(t, PredF.PVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h3,h4], leg_act, 'Representation Model' , 'Pos','Delayed Vel');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    
    % arrange for plotErrorBarGain
    g_PVtau_lateAdapt=[K_PVtau_lateAdapt(g) Bt_PVtau_lateAdapt(g)];
    gint_PVtau_lateAdapt=[Kint_PVtau_lateAdapt(g,:)' Btint_PVtau_lateAdapt(g,:)'];
    gs_PVtau_lateAdapt=[Ks_PVtau_lateAdapt{g} Bts_PVtau_lateAdapt{g}];
    c_PVt=[cP;cVt];
    legGain={'Pos','Delayed Vel'};
    
    figure('Position',fs_bar3)
    hold on
    plotErrorBarGain(g_PVtau_lateAdapt,gint_PVtau_lateAdapt,gs_PVtau_lateAdapt,c_PVt,legGain)
    
    % VA
    PredF.VA.BV.all=VFactor*B_VA_lateAdapt(g)*VelY.all;
    PredF.VA.MA.all=AFactor*M_VA_lateAdapt(g)*AccY.all;
    PredF.VA.Tot.all=PredF.VA.BV.all + PredF.VA.MA.all;
    PredF.VA.BV.m=nanmean(PredF.VA.BV.all); % mean of all subjects
    PredF.VA.BV.ste=nanstd(PredF.VA.BV.all)/sqrt(trForSTE);
    PredF.VA.BV.ci=NonNanDataCI(PredF.VA.BV.all);
    PredF.VA.MA.m=nanmean(PredF.VA.MA.all); % mean of all subjects
    PredF.VA.MA.ste=nanstd(PredF.VA.MA.all)/sqrt(trForSTE);
    PredF.VA.MA.ci=NonNanDataCI(PredF.VA.MA.all);
    PredF.VA.Tot.m=nanmean(PredF.VA.Tot.all); % mean of all subjects
    PredF.VA.Tot.ste=nanstd(PredF.VA.Tot.all)/sqrt(trForSTE);
    PredF.VA.Tot.ci=NonNanDataCI(PredF.VA.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p4=plot(t, PredF.VA.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    p5=plot(t, PredF.VA.MA.m,':', 'color' ,cA , 'linewidth' ,lw);
    p1=plot(t, PredF.VA.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,PredF.VA.Tot.ci(2,:),PredF.VA.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,PredF.VA.BV.ci(2,:),PredF.VA.BV.ci(1,:),1-((1-cV)*scJb));
    jbfill(t,PredF.VA.MA.ci(2,:),PredF.VA.MA.ci(1,:),1-((1-cA)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h4=plot(t, PredF.VA.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    h5=plot(t, PredF.VA.MA.m,':', 'color' ,cA , 'linewidth' ,lw);
    h1=plot(t, PredF.VA.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h4,h5], leg_act, 'Representation Model' ,'Vel','Acc');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    
    % arrange for plotErrorBarGain
    g_VA_lateAdapt=[B_VA_lateAdapt(g) M_VA_lateAdapt(g)];
    gint_VA_lateAdapt=[Bint_VA_lateAdapt(g,:)' Mint_VA_lateAdapt(g,:)'];
    gs_VA_lateAdapt=[Bs_VA_lateAdapt{g} Ms_VA_lateAdapt{g}];
    c_VA=[cV;cA];
    legGain={'Vel','Acc'};
    
    figure('Position',fs_bar3)
    hold on
    plotErrorBarGain(g_VA_lateAdapt,gint_VA_lateAdapt,gs_VA_lateAdapt,c_VA,legGain)
    
    % VVtau
    PredF.VVtau.BV.all=VFactor*B_VVtau_lateAdapt(g)*VelY.all;
    PredF.VVtau.BtVtau.all=VFactor*Bt_VVtau_lateAdapt(g)*VeltauY.all;
    PredF.VVtau.Tot.all=PredF.VVtau.BV.all + PredF.VVtau.BtVtau.all;
    PredF.VVtau.BV.m=nanmean(PredF.VVtau.BV.all); % mean of all subjects
    PredF.VVtau.BV.ste=nanstd(PredF.VVtau.BV.all)/sqrt(trForSTE);
    PredF.VVtau.BV.ci=NonNanDataCI(PredF.VVtau.BV.all);
    PredF.VVtau.BtVtau.m=nanmean(PredF.VVtau.BtVtau.all); % mean of all subjects
    PredF.VVtau.BtVtau.ste=nanstd(PredF.VVtau.BtVtau.all)/sqrt(trForSTE);
    PredF.VVtau.BtVtau.ci=NonNanDataCI(PredF.VVtau.BtVtau.all);
    PredF.VVtau.Tot.m=nanmean(PredF.VVtau.Tot.all); % mean of all subjects
    PredF.VVtau.Tot.ste=nanstd(PredF.VVtau.Tot.all)/sqrt(trForSTE);
    PredF.VVtau.Tot.ci=NonNanDataCI(PredF.VVtau.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p4=plot(t, PredF.VVtau.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    p5=plot(t, PredF.VVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    p1=plot(t, PredF.VVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,PredF.VVtau.Tot.ci(2,:),PredF.VVtau.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,PredF.VVtau.BV.ci(2,:),PredF.VVtau.BV.ci(1,:),1-((1-cV)*scJb));
    jbfill(t,PredF.VVtau.BtVtau.ci(2,:),PredF.VVtau.BtVtau.ci(1,:),1-((1-cVt)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h4=plot(t, PredF.VVtau.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    h5=plot(t, PredF.VVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    h1=plot(t, PredF.VVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h4,h5], leg_act, 'Representation Model' ,'Vel','Delayed Vel');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    
    % arrange for plotErrorBarGain
    g_VVtau_lateAdapt=[B_VVtau_lateAdapt(g) Bt_VVtau_lateAdapt(g)];
    gint_VVtau_lateAdapt=[Bint_VVtau_lateAdapt(g,:)' Btint_VVtau_lateAdapt(g,:)'];
    gs_VVtau_lateAdapt=[Bs_VVtau_lateAdapt{g} Bts_VVtau_lateAdapt{g}];
    c_VVtau=[cV;cVt];
    legGain={'Vel','Delayed Vel'};
    
    figure('Position',fs_bar3)
    hold on
    plotErrorBarGain(g_VVtau_lateAdapt,gint_VVtau_lateAdapt,gs_VVtau_lateAdapt,c_VVtau,legGain)
    
    % PVA
    PredF.PVA.KP.all=PFactor*K_PVA_lateAdapt(g)*PosY.all;
    PredF.PVA.BV.all=VFactor*B_PVA_lateAdapt(g)*VelY.all;
    PredF.PVA.MA.all=AFactor*M_PVA_lateAdapt(g)*AccY.all;
    PredF.PVA.Tot.all=PredF.PVA.KP.all + PredF.PVA.BV.all + PredF.PVA.MA.all;
    PredF.PVA.KP.m=nanmean(PredF.PVA.KP.all); % mean of all subjects
    PredF.PVA.KP.ste=nanstd(PredF.PVA.KP.all)/sqrt(trForSTE);
    PredF.PVA.KP.ci=NonNanDataCI(PredF.PVA.KP.all);
    PredF.PVA.BV.m=nanmean(PredF.PVA.BV.all); % mean of all subjects
    PredF.PVA.BV.ste=nanstd(PredF.PVA.BV.all)/sqrt(trForSTE);
    PredF.PVA.BV.ci=NonNanDataCI(PredF.PVA.BV.all);
    PredF.PVA.MA.m=nanmean(PredF.PVA.MA.all); % mean of all subjects
    PredF.PVA.MA.ste=nanstd(PredF.PVA.MA.all)/sqrt(trForSTE);
    PredF.PVA.MA.ci=NonNanDataCI(PredF.PVA.MA.all);
    PredF.PVA.Tot.m=nanmean(PredF.PVA.Tot.all); % mean of all subjects
    PredF.PVA.Tot.ste=nanstd(PredF.PVA.Tot.all)/sqrt(trForSTE);
    PredF.PVA.Tot.ci=NonNanDataCI(PredF.PVA.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p3=plot(t, PredF.PVA.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    p4=plot(t, PredF.PVA.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    p5=plot(t, PredF.PVA.MA.m,':', 'color' ,cA , 'linewidth' ,lw);
    p1=plot(t, PredF.PVA.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,PredF.PVA.Tot.ci(2,:),PredF.PVA.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,PredF.PVA.KP.ci(2,:),PredF.PVA.KP.ci(1,:),1-((1-cP)*scJb));
    jbfill(t,PredF.PVA.BV.ci(2,:),PredF.PVA.BV.ci(1,:),1-((1-cV)*scJb));
    jbfill(t,PredF.PVA.MA.ci(2,:),PredF.PVA.MA.ci(1,:),1-((1-cA)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h3=plot(t, PredF.PVA.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    h4=plot(t, PredF.PVA.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    h5=plot(t, PredF.PVA.MA.m,':', 'color' ,cA , 'linewidth' ,lw);
    h1=plot(t, PredF.PVA.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h3,h4,h5], leg_act, 'Representation Model' ,'Pos','Vel','Acc');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    
    % arrange for plotErrorBarGain
    g_PVA_lateAdapt=[K_PVA_lateAdapt(g) B_PVA_lateAdapt(g) M_PVA_lateAdapt(g)];
    gint_PVA_lateAdapt=[Kint_PVA_lateAdapt(g,:)' Bint_PVA_lateAdapt(g,:)' Mint_PVA_lateAdapt(g,:)'];
    gs_PVA_lateAdapt=[Ks_PVA_lateAdapt{g} Bs_PVA_lateAdapt{g} Ms_PVA_lateAdapt{g}];
    c_PVA=[cP;cV;cA];
    legGain={'Pos','Vel','Acc'};
    
    figure('Position',fs_bar3)
    hold on
    plotErrorBarGain(g_PVA_lateAdapt,gint_PVA_lateAdapt,gs_PVA_lateAdapt,c_PVA,legGain)
    
    % PVVtau
    PredF.PVVtau.KP.all=PFactor*K_PVVtau_lateAdapt(g)*PosY.all;
    PredF.PVVtau.BV.all=VFactor*B_PVVtau_lateAdapt(g)*VelY.all;
    PredF.PVVtau.BtVtau.all=VFactor*Bt_PVVtau_lateAdapt(g)*VeltauY.all;
    PredF.PVVtau.Tot.all=PredF.PVVtau.KP.all + PredF.PVVtau.BV.all + PredF.PVVtau.BtVtau.all;
    PredF.PVVtau.KP.m=nanmean(PredF.PVVtau.KP.all); % mean of all subjects
    PredF.PVVtau.KP.ste=nanstd(PredF.PVVtau.KP.all)/sqrt(trForSTE);
    PredF.PVVtau.KP.ci=NonNanDataCI(PredF.PVVtau.KP.all);
    PredF.PVVtau.BV.m=nanmean(PredF.PVVtau.BV.all); % mean of all subjects
    PredF.PVVtau.BV.ste=nanstd(PredF.PVVtau.BV.all)/sqrt(trForSTE);
    PredF.PVVtau.BV.ci=NonNanDataCI(PredF.PVVtau.BV.all);
    PredF.PVVtau.BtVtau.m=nanmean(PredF.PVVtau.BtVtau.all); % mean of all subjects
    PredF.PVVtau.BtVtau.ste=nanstd(PredF.PVVtau.BtVtau.all)/sqrt(trForSTE);
    PredF.PVVtau.BtVtau.ci=NonNanDataCI(PredF.PVVtau.BtVtau.all);
    PredF.PVVtau.Tot.m=nanmean(PredF.PVVtau.Tot.all); % mean of all subjects
    PredF.PVVtau.Tot.ste=nanstd(PredF.PVVtau.Tot.all)/sqrt(trForSTE);
    PredF.PVVtau.Tot.ci=NonNanDataCI(PredF.PVVtau.Tot.all);
    
    figure('Position',fs)
    hold on
    p2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    p3=plot(t, PredF.PVVtau.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    p4=plot(t, PredF.PVVtau.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    p5=plot(t, PredF.PVVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    p1=plot(t, PredF.PVVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    jbfill(t,PredF.PVVtau.Tot.ci(2,:),PredF.PVVtau.Tot.ci(1,:),1-((1-c_rep)*scJb));
    jbfill(t,ActF.ci(2,:),ActF.ci(1,:),1-((1-c_act)*scJb));
    jbfill(t,PredF.PVVtau.KP.ci(2,:),PredF.PVVtau.KP.ci(1,:),1-((1-cP)*scJb));
    jbfill(t,PredF.PVVtau.BV.ci(2,:),PredF.PVVtau.BV.ci(1,:),1-((1-cV)*scJb));
    jbfill(t,PredF.PVVtau.BtVtau.ci(2,:),PredF.PVVtau.BtVtau.ci(1,:),1-((1-cVt)*scJb));
    hold on
    h2=plot(t, ActF.m,'--','color',c_act, 'linewidth' ,lw); % data
    h3=plot(t, PredF.PVVtau.KP.m,':', 'color' ,cP , 'linewidth' ,lw);
    h4=plot(t, PredF.PVVtau.BV.m,':', 'color' ,cV , 'linewidth' ,lw);
    h5=plot(t, PredF.PVVtau.BtVtau.m,':', 'color' ,cVt , 'linewidth' ,lw);
    h1=plot(t, PredF.PVVtau.Tot.m, 'color' ,c_rep  , 'linewidth' ,lwRep);
    line([0  UniteTrLen(g)*st*s2ms] , [0  0] , 'Color' , [0.4  0.4  0.4] ,'linestyle' , ':'  , 'linewidth' , 2);
    ylim(yLim);
    xlim(xLim);
    xlabel('Time [ms]','fontsize' , 20, 'fontweight' , 'b');
    ylabel('Force [N]','fontsize' , 20, 'fontweight' , 'b');
    l=legend([h2,h1,h3,h4,h5], leg_act, 'Representation Model' ,'Pos','Vel','Delayed Vel');
    set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
    legend boxoff
    set(gca,'box','off', 'FontSize', afs,'fontweight','b','xtick',xtick,'xticklabel',xticklabel,'ytick',ytick)
    
    % arrange for plotErrorBarGain
    g_PVVtau_lateAdapt=[K_PVVtau_lateAdapt(g) B_PVVtau_lateAdapt(g) Bt_PVVtau_lateAdapt(g)];
    gint_PVVtau_lateAdapt=[Kint_PVVtau_lateAdapt(g,:)' Bint_PVVtau_lateAdapt(g,:)' Btint_PVVtau_lateAdapt(g,:)'];
    gs_PVVtau_lateAdapt=[Ks_PVVtau_lateAdapt{g} Bs_PVVtau_lateAdapt{g} Bts_PVVtau_lateAdapt{g}];
    c_PVVtau=[cP;cV;cVt];
    legGain={'Pos','Vel','Delayed Vel'};
    
    figure('Position',fs_bar3)
    hold on
    plotErrorBarGain(g_PVVtau_lateAdapt,gint_PVVtau_lateAdapt,gs_PVVtau_lateAdapt,c_PVVtau,legGain)
    
    if g==1
        
    elseif g==2
        % save the gains of the primitives in this group from the PVA and
        % PVVtau models - for the predictions of slow to fast
        % generalization
        save('g_D70_lateAdapt_PVA_PVVtau','g_PVA_lateAdapt','g_PVVtau_lateAdapt');
    elseif g==3
        
    else
        save('gs_D70SF_lateAdapt_PVA_PVVtau_VA_VVtau','gs_PVA_lateAdapt','gs_PVVtau_lateAdapt','gs_VA_lateAdapt','gs_VVtau_lateAdapt');
        maxVelY=max(nanmean(VelY.all));
        maxActF=max(ActF.m);
        ActF_VelY_factor_D70SF_LA=maxActF/maxVelY; 
        save('ActF_VelY_factor_D70SF_LA','ActF_VelY_factor_D70SF_LA');
    end
    
end
