function [p,statobs,statperm] = permtest(stat,y,y2,nrand,diffflag)

if nargin < 5
   diffflag = 0;
end

if nargin < 4
   nrand = 1000; 
end

if isvector(y)
    y = vert(y);
end

if isvector(y2)
    y2 = vert(y2);
end


if diffflag && isa(stat,'function_handle')
   stat = @(x1,x2)(stat(x1)-stat(x2));
end

if ischar(stat)
   switch stat
       case {'anova1'}     
           stat = @(y,grp)anova1(y,grp,'off');
   end
end

statobs = stat(y,y2);

if diffflag == 1
   n1 = length(y); n2 = length(y2);
   y = cat(1,vert(y),vert(y2));
   y2 = Make_designVect([n1 n2])';
end

for i = 1:nrand
    perm = randperm(size(y,1));
    permy2 = y2(perm,:);
    if diffflag == 1
        statperm(i) = stat(y(permy2==1,:),y(permy2==2,:));
    else
        statperm(i) = stat(y,permy2);
    end
end

p = value_prctile(statperm,statobs);
if p > 0.5
    p = 1-p;
end
p = p*2; %two tailed test