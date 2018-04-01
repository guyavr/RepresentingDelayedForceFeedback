function plotStatForceRegAdapt_ND_D70_D100(mean,err)

c = get(gcf,'DefaultAxesColorOrder');
% close

cnd=c(1,:)*1.2;
cd70=c(3,:)*1.05;
cd100=c(7,:)*1.05;

lw=5;
eblw=2;
bw=0.25; % bar width
pw=bw/2; % patch half widht
bd=0.3; % bars' distances of the same stage
afs=26;
tickfs=24;
astfs=28; % fontsize of asterisk
y_lim=[0 1.15];
x_lim=[0.5 2.5];

scrsz = [1,1,1280,800];
figure('Position',[scrsz(3)/4 scrsz(4)*1/8 scrsz(3)*2.1/4 scrsz(4)*2.5/4])
hold on
b(1)=bar(1-bd,mean(1,1),bw);
b(2)=bar(2-bd,mean(1,2),bw);
b(3)=bar(1,mean(2,1),bw);
b(4)=bar(2,mean(2,2),bw);
b(5)=bar(1+bd,mean(3,1),bw);
b(6)=bar(2+bd,mean(3,2),bw);
box off
set(b(1),'facecolor','none','edgecolor',cnd,'linewidth',lw);
set(b(2),'facecolor',cnd,'edgecolor',cnd,'linewidth',lw);
set(b(3),'facecolor','none','edgecolor',cd70,'linewidth',lw);
set(b(4),'facecolor',cd70,'edgecolor',cd70,'linewidth',lw);
set(b(5),'facecolor','none','edgecolor',cd100,'linewidth',lw);
set(b(6),'facecolor',cd100,'edgecolor',cd100,'linewidth',lw);
l1=patch([1-bd-pw 1-bd+pw 1-bd+pw 1-bd-pw],[0 0 mean(1,1) mean(1,1)],[0 0 0 0],'linewidth',lw, 'edgecolor',cnd,'markerfacecolor','flat','facecolor','w');
hatch(l1,45,cnd,'-',18,3);
l2=patch([1-pw 1+pw 1+pw 1-pw],[0 0 mean(2,1) mean(2,1)],[0 0 0 0],'linewidth',lw, 'edgecolor',cd70,'markerfacecolor','flat','facecolor','w');
hatch(l2,45,cd70,'-',18,3);
l3=patch([1+bd-pw 1+bd+pw 1+bd+pw 1+bd-pw],[0 0 mean(3,1) mean(3,1)],[0 0 0 0],'linewidth',lw, 'edgecolor',cd100,'markerfacecolor','flat','facecolor','w');
hatch(l3,45,cd100,'-',18,3);

% xlabel('Stage','fontsize',26,'fontweight' , 'b');
% set(gca,'layer','bottom','XTick',1:2,'XTickLabel',{'EA','LA'},'ytick',[0.2:0.2:1],'fontsize',tickfs,'fontweight','b','box','off');
set(gca,'layer','bottom','XTick',1:2,'XTickLabel',{'     Early\newlineAdaptation','     Late\newlineAdaptation'},'ytick',[0.2:0.2:1],'fontsize',tickfs,'fontweight','b','box','off');
ylabel('Mean Adaptation Coeff.','fontsize',afs-2,'fontweight' , 'b');
% l=legend ([b(1),b(5)],'Non-delayed','Delayed','location','southwest');
% set(l, 'Box', 'off','fontsize',21,'fontweight' , 'b')

xlim(x_lim);
ylim(y_lim);

errup=mean+err;
e(1)=line([1 1]-bd,[mean(1,1) errup(1,1)],'color',cnd,'linewidth',lw);
e(2)=line([2 2]-bd,[mean(1,2) errup(1,2)],'color',cnd,'linewidth',lw);
e(3)=line([1 1],[mean(2,1) errup(2,1)],'color',cd70,'linewidth',lw);
e(4)=line([2 2],[mean(2,2) errup(2,2)],'color',cd70,'linewidth',lw);
e(5)=line([1 1]+bd,[mean(3,1) errup(3,1)],'color',cd100,'linewidth',lw);
e(6)=line([2 2]+bd,[mean(3,2) errup(3,2)],'color',cd100,'linewidth',lw);
% hold on
ebwid=0.06;
line([1-bd-ebwid  1-bd+ebwid] , [errup(1,1) errup(1,1)] , 'color' , cnd , 'linewidth' , lw);
line([2-bd-ebwid  2-bd+ebwid] , [errup(1,2) errup(1,2)] , 'color' , cnd , 'linewidth' , lw);
line([1-ebwid  1+ebwid] , [errup(2,1) errup(2,1)] , 'color' , cd70 , 'linewidth' , lw);
line([2-ebwid  2+ebwid] , [errup(2,2) errup(2,2)] , 'color' , cd70 , 'linewidth' , lw);
line([1+bd-ebwid  1+bd+ebwid] , [errup(3,1) errup(3,1)] , 'color' , cd100 , 'linewidth' , lw);
line([2+bd-ebwid  2+bd+ebwid] , [errup(3,2) errup(3,2)] , 'color' , cd100 , 'linewidth' , lw);

% % significance asterisk - 11Jun16 with the significance between D100 and the other
% % groups during LA
% lineLoc=[0.8,0.86,0.98,1.04,1.1];
% astLoc=lineLoc+0.01;
% line([2 2+bd],[lineLoc(1) lineLoc(1)],'color','k','linewidth',eblw);
% text(2+bd/2,astLoc(1),'**','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
% line([2-bd 2+bd],[lineLoc(2) lineLoc(2)],'color','k','linewidth',eblw);
% text(2,astLoc(2),'**','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
% line([1+bd 2+bd],[lineLoc(3) lineLoc(3)],'color','k','linewidth',eblw);
% text(1.5+bd,astLoc(3),'***','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
% line([1 2],[lineLoc(4) lineLoc(4)],'color','k','linewidth',eblw);
% text(1.5,astLoc(4),'***','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
% line([1-bd 2-bd],[lineLoc(5) lineLoc(5)],'color','k','linewidth',eblw);
% text(1.5-bd,astLoc(5),'***','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');

% significance asterisk - 11Jun16 - no interaction effect
lineLoc=0.2+[0.8,0.87,0.94,0.67,0.59];
astLoc=lineLoc+0.01;
line([1+bd 2+bd],[lineLoc(1) lineLoc(1)],'color','k','linewidth',eblw);
text(1.5+bd,astLoc(1),'***','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
line([1 2],[lineLoc(2) lineLoc(2)],'color','k','linewidth',eblw);
text(1.5,astLoc(2),'***','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
line([1-bd 2-bd],[lineLoc(3) lineLoc(3)],'color','k','linewidth',eblw);
text(1.5-bd,astLoc(3),'***','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
line([2-bd 2+bd],[lineLoc(4) lineLoc(4)],'color','k','linewidth',eblw);
text(2,astLoc(4),'**','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
line([2 2+bd],[lineLoc(5) lineLoc(5)],'color','k','linewidth',eblw);
text(2+bd/2,astLoc(5),'**','color','k','fontsize',astfs,'fontweight','b','horizontalalignment','center');
end

