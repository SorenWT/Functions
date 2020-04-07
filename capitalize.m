function capcell = capitalize(cellin)

for i = 1:length(cellin)
   capcell{i} = [upper(cellin{i}(1)) cellin{i}(2:end)];
end