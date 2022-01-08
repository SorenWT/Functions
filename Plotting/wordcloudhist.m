function wordcloudhist(data,reorder)

tmp = wordCloudCounts(data);

if nargin < 2
   reorder = tmp.Word;
end

[~,m1] = match_str(tmp.Word,reorder);
tmp = tmp(m1,:);

bar(tmp.Count)
set(gca,'XTick',1:length(tmp.Count),'XTickLabel',tmp.Word)
FixAxes(gca,16)
xtickangle(90)
ylabel('Count')
