function [R,P,c,t] = imagesc_stars(R,P,labels,fsize,dobelowdiag)
% plotCorrelationMatrix
% Plot correlation matrix
%
% plotCorrelationMatrix(X,labels) creates a plot of the matrix of Pearson
% correlation coefficients between the variables in X.
% plotCorrelationMatrix displays the correlations between each two
% variables as a colored cell in a table, ranging from red for positive
% correlations to white for zero correlations, to blue for negative
% correlations. The correlation coefficient is displayed in each cell with
% asterisks indicating the p value. Only the below-diagonal correlations
% are displayed for simplicity.
%
% [R,P] = plotCorrelationMatrix(X,labels) also returns the matrix of
% correlation coefficients and the matrix of p-values for testing the
% hypothesis that there is no relationship between the variables (null
% hypothesis).

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

cla reset % clear and reset figure axes

nVar = size(R,2);

% if labels are not provided, use numbers
if nargin==1
    labels = num2cell(1:nVar);
end

if nargin < 4
   fsize = 18; %swt edit for changing font size 
end

if nargin < 5
   dobelowdiag = 1; 
end

if size(R,1) ~= size(R,2)
   dobelowdiag = 0; 
end

% compute correlations
%[R,P] = corrcoef(X,'rows','pairwise');

% hide above-diagonal values by making r values above the diagonal 0 and
% p values above the diagonal 1
[rowIndex,colIndex] = find(~isnan(R));
if dobelowdiag
R(colIndex>rowIndex) = 0;
P(colIndex>rowIndex) = 1;
end

% create colormap
blueWhiteRed = [createColorGradient([0 0 1],[1 1 1],32);...
    createColorGradient([1 1 1],[1 0 0],32)];

% create image
Rmask = zeros(size(R));
Rmask(~isnan(R)) = R(~isnan(R));
imagesc(Rmask)
colormap(gca,blueWhiteRed)
c = colorbar('eastoutside');
clim = get(gca,'CLim');
maxlim = max(abs(clim));
caxis([-maxlim maxlim])
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','none')
set(gca,'XTick',1:size(R,2))
if dobelowdiag
set(gca,'XTickLabels',labels)
end
nChar = cellfun(@length,labels);
if any(nChar>1)
    set(gca,'XTickLabelRotation',45)
end
set(gca,'YTick',1:size(R,1))
if dobelowdiag
set(gca,'YTicklabels',labels)
end
set(gca,'TickLength',[0 0])
axis square; box off

% add text displaying r-values with asterisks for significance
for row = 1 : size(R,1)
    for col = 1 : size(R,2)
        if ~dobelowdiag || row>col % below diagonal
            %str = num2str(R(row,col),2);
            if P(row,col)<0.001
                str = ['***'];
            elseif P(row,col)<0.01
                str = ['**'];
            elseif P(row,col)<0.05
                str = ['*'];
            elseif P(row,col)<0.1
                str = ['\diamondsuit'];
            else
                str = [' '];
            end
            text(col,row,str,'HorizontalAlignment','center','FontSize',fsize);
        end
    end
end

% add legend for significance
if dobelowdiag
t = text(nVar*0.8,1,{'\diamondsuit   p<0.1','*   p<0.05','**  p<0.01','*** p<0.001'},...
    'HorizontalAlignment','Left');
else
   t = text(nVar+0.5,0.5,{'\diamondsuit   p<0.1','*   p<0.05','**  p<0.01','*** p<0.001'},...
    'HorizontalAlignment','Left');
end
uistack(t,'top')

end

function colorGradient = createColorGradient(startColor,endColor,nColors)

colorGradient = zeros(nColors,3);
for i = 1 : 3
    if endColor(i)~=startColor(i)
        colorGradient(:,i) = startColor(i) : (endColor(i)-startColor(i))/(nColors-1) : endColor(i);
    elseif endColor(i)==startColor(i)
        colorGradient(:,i) = startColor(i);
    end
end

end