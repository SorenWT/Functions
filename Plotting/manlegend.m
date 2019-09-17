function manlegend(labels,colors)

for c = 1:size(colors,1)
    patch([NaN NaN],[NaN NaN],colors(c,:))
end

legend(labels)