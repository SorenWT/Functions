function counts = wordcloudhist(data,reorder)

counts = wordCloudCounts(data);

if nargin < 2
   reorder = counts.Word;
end

[~,m1] = match_str(counts.Word,reorder);
counts = counts(m1,:);

bar(counts.Count)
set(gca,'XTick',1:length(counts.Count),'XTickLabel',counts.Word)
FixAxes(gca,16)
xtickangle(90)
ylabel('Count')
