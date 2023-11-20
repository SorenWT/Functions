function [x,pdf] = probplot(data)

[y,x] = ksdensity(data);
hold on
plot(x,y,'LineWidth',2)