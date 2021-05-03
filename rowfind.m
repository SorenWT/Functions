function indxout = rowfind(arrin,n,dir)

if dir == 1
   arrin = arrin';
end

for i = 1:size(arrin,1)
    tmp = find(arrin(i,:),n); 
    if ~isempty(tmp)
        indxout(i) = tmp;
    else
        indxout(i) = NaN;
    end
end