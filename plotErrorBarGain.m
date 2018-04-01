function plotErrorBarGain(g,gint,gs,c,lg)
% Plot the bars of the primitives gain (2 is for 2 primitives)
% inputs:
% g- a vector with the value of each gain
% gint- matrix- each column is for each primitive, 1st row is lower ci, 2nd
% row is the upper ci
% c- colors
% lg- cell array of strings for the legend

afs=19;
lfs=18;

ebwid=0.1; % horizontal line errorbar width
elw=3; % errorbar linewidth
lw=3;
ms=7;

eb_shift=.2;

ybar_lim=[-0.32 1.18];
xbar_lim=[0.5 3.5];

Ns=size(gs,1); % number of subjects
Np=length(g); % number of primitives
for i=1:Np
    b(i)=bar(i,g(i));
    set(b(i),'facecolor',c(i,:),'edgecolor',c(i,:),'linewidth',lw);
    e(i)=errorbar(i-eb_shift,g(i),g(i)-gint(1,i),gint(2,i)-g(i),'Color',c(i,:),'linewidth',elw);
    if g(i)>0
        line([i i]-eb_shift , [g(i) gint(1,i)],'Color','w', 'linewidth' , elw);
        line([i-ebwid  i+ebwid]-eb_shift , [gint(1,i) gint(1,i)],'Color','w', 'linewidth' , elw);
        line([i-ebwid  i+ebwid]-eb_shift , [gint(2,i) gint(2,i)],'Color' ,c(i,:) , 'linewidth' , elw);
        for s=1:Ns
            if (gs(s,i)>0 && gs(s,i)<g(i))
                plot(i+eb_shift,gs(s,i),'o','markerfacecolor','w','markeredgecolor','k','markersize',ms)
            else
                plot(i+eb_shift,gs(s,i),'o','markerfacecolor',c(i,:),'markeredgecolor','k','markersize',ms)
            end
        end
    else
        line([i i]-eb_shift , [g(i) gint(2,i)],'Color','w', 'linewidth' , elw);
        line([i-ebwid  i+ebwid]-eb_shift , [gint(1,i) gint(1,i)],'Color',c(i,:) , 'linewidth' , elw);
        line([i-ebwid  i+ebwid]-eb_shift , [gint(2,i) gint(2,i)],'Color' ,'w' , 'linewidth' , elw);
        for s=1:Ns
            if (gs(s,i)<0 && gs(s,i)>g(i))
                plot(i+eb_shift,gs(s,i),'o','markerfacecolor','w','markeredgecolor','k','markersize',ms)
            else
                plot(i+eb_shift,gs(s,i),'o','markerfacecolor',c(i,:),'markeredgecolor','k','markersize',ms)
            end
        end
%         % only for Exp2Vel - PVA model during late adaptation
%         line([i i]-eb_shift , [g(i) gint(2,i)],'Color',c(i,:), 'linewidth' , elw);
%         line([i-ebwid  i+ebwid] , [gint(1,i) gint(1,i)],'Color',c(i,:) , 'linewidth' , elw);
%         line([i-ebwid  i+ebwid] , [gint(2,i) gint(2,i)],'Color' ,c(i,:) , 'linewidth' , elw);
    end
    
    
end
set(gca,'layer','bottom','XTick',1:Np,'XTickLabel',[],'ytick',[-0.2:0.2:1],'fontsize',afs,'fontweight','b','box','off');
box off
ylabel('Gain','fontsize',20,'fontweight' , 'b');
xlim(xbar_lim);
ylim(ybar_lim);
l=legend(b,lg);
set(l, 'fontweight' , 'b','fontsize' , lfs,'location','northwest');
legend boxoff

end

