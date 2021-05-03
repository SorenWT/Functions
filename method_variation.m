function mvar = method_variation(stats,pfield)

count = 1;
tmpstats = stats;
while 1
    fields{count} = fieldnames(tmpstats);
    try
        getfield_nest(tmpstats,pfield);
        fields(end) = [];
        count = count-1;
        break;
    catch
                tmpstats = tmpstats.(fields{count}{1});
        count = count+1;
    end
%     
%     if isfield_nest(tmpstats,pfield) 
%         fields(end) = [];
%         count = count-1;
%         break;
%     else
%         tmpstats = tmpstats.(fields{count}{1});
%         count = count+1;
%     end
end

if count == 0
    error('You must supply a nested stats structure with more than one version of the analysis')
end

fields = fieldnames_recurse(stats);
fields = cell_unpack(fields);

fields = cellfun(@(d)tokenize(d,'.'),fields,'UniformOutput',false);

for i = 1:length(fields)
    fields(i) = join(fields{i}(1:count),'.');
end
fields = unique(fields);

for i = 1:length(fields)
    p{i} = getfield_nest(stats,[fields{i} '.' pfield]);
    p{i} = vert(p{i}); % in case there are multiple p values
end

p = cat(2,p{:});

mvar.pall = array2table(p,'VariableNames',fields);
mvar.pvar = nanstd(p,[],2);
mvar.sigindx05 = p < 0.05; mvar.sigindx10 = p < 0.1;
mvar.prcsig05 = nanmean(p < 0.05,2); mvar.prcsig10 = nanmean(p < 0.1,2);


