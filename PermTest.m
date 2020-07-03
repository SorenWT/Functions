function p = permtest(stat,y,group,nrand,diffflag)

if nargin < 5
   diffflag = 0;
end

if nargin < 4
   nrand = 1000; 
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

statobs = stat(y,group);

if diffflag == 1
   n1 = length(y); n2 = length(group);
   y = cat(1,vert(y),vert(group));
   group = Make_designVect([n1 n2])';
end

for i = 1:nrand
    perm = randperm(length(y));
    permgroup = group(perm);
    if diffflag == 1
        statperm(i) = stat(y(permgroup==1),y(permgroup==2));
    else
        statperm(i) = stat(y,permgroup);
    end
end

p = value_prctile(statperm,statobs);
if p > 0.5
    p = 1-p;
end
p = p*2; %two tailed test