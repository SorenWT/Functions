function rsaplot(rsadata,datatitles,maintitle,p,pindx)

if ~exist('p','var')
    p = panel('no-manage-font');
    pindx = {};
end

p(pindx{:}).pack('h',{0.3 0.4 0.3});

if ~exist('maintitle','var') || isempty(maintitle)
    maintitle = 'Representational similarity';
end

if ~exist('datatitles','var') || isempty(datatitles)
    datatitles = {'Data 1','Data 2'};
end

for i = 1:2
    if strcmpi(rsadata.type{i},'pearson') || strcmpi(rsadata.type{i},'spearman')
        rsadata.corrmats{i} = 1-rsadata.corrmats{i};
        rsadata.corrvals{i} = 1-rsadata.corrvals{i};
    end
end

anynan = find(isnan(rsadata.corrmats{1}(:,1)) | isnan(rsadata.corrmats{2}(:,1)));

for i = 1:2
        rsadata.corrmats{i}(:,anynan) = [];
        rsadata.corrmats{i}(anynan,:) = [];
end


p(pindx{:},1).select()
imagesc(rsadata.corrmats{1})
title([datatitles{1} ' distance matrix'],'FontSize',18)
cbar = colorbar('westoutside');
switch rsadata.type{1}
    case 'pearson'
        cbar.Label.String = "1 - Pearson's r";
    case 'spearman'
        cbar.Label.String = "1 - Spearman's \rho";
    case 'eucdist'
        cbar.Label.String = "log(euclidean distance)";
end


cbar.FontSize = 14;
set(gca,'XLim',[0.5 size(rsadata.corrmats{1},1)+0.5],'YLim',...
    [0.5 size(rsadata.corrmats{1},1)+0.5],'YDir','reverse')
axis square

p(pindx{:},3).select()
imagesc(rsadata.corrmats{2})
title([datatitles{2} ' distance matrix'],'FontSize',18)
cbar = colorbar;

switch rsadata.type{2}
    case 'pearson'
        cbar.Label.String = "1 - Pearson's r";
    case 'spearman'
        cbar.Label.String = "1 - Spearman's \rho";
    case 'eucdist'
        cbar.Label.String = "Euclidean distance";
end
cbar.FontSize = 14;
set(gca,'XLim',[0.5 size(rsadata.corrmats{2},1)+0.5],'YLim',...
    [0.5 size(rsadata.corrmats{2},1)+0.5],'YDir','reverse')
axis square


p(pindx{:},2).select()



nicecorrplot(rsadata.corrvals{1},rsadata.corrvals{2},{[datatitles{1} ' distance'],[datatitles{2} ' distance']},...
    'externalp',rsadata.pperm,'type','pearson');
title(maintitle,'FontSize',22)

p(pindx{:},3).marginleft = 8;

if isempty(pindx)
    p.margintop = 10;
    p.de.marginleft = 30;
    p.marginleft = 30;
    p.marginright = 25;
end