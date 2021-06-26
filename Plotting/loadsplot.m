function loadsplot(loads,ci,clr,plottype)
% rows of CI should be lower / upper

if nargin < 3
    l = lines;
    clr = l(1,:);
end

if nargin < 4
    plottype = 'bar';
end

if isvector(ci)
    ci = [horz(loads)-horz(ci); horz(loads)+horz(ci)];
end

if ~all(size(ci)==2) && size(ci,2)==2
    ci = ci';
end

switch plottype
    case 'bar'
        bar(1:length(loads),loads,'facecolor',palecol(clr));
        hold on
        bar(1:length(loads),loads.*((ci(1,:)>0 & ci(2,:)>0) | (ci(1,:)<0 & ci(2,:)<0)),'facecolor',clr)
        errorbar(1:length(loads),loads,loads-ci(1,:),ci(2,:)-loads,'linestyle','none','linewidth',2,'color','k')
    case 'barh'
          barh(1:length(loads),loads,'facecolor',palecol(clr));
        hold on
        barh(1:length(loads),loads.*((ci(1,:)>0 & ci(2,:)>0) | (ci(1,:)<0 & ci(2,:)<0)),'facecolor',clr)
        errorbar(loads,1:length(loads),loads-ci(1,:),ci(2,:)-loads,'horizontal','linestyle','none','linewidth',2,'color','k')
      
end

