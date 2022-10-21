function [b,err] = errorbar_bar(x,y,err)

if isempty(x)
   x = 1:length(y); 
end

% quick function for doing bar plots with errorbars

b = bar(y);
hold on

for i = 1:length(b)
    errorbar(b(i).XEndPoints,y(:,i),err(i,:),'LineStyle','none','LineWidth',3,'color','k')
end

FixAxes(gca,18)

