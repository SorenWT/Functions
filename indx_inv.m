function inv = indx_inv(indx)

for i = 1:length(indx)
    inv(i) = find(indx==i); 
end