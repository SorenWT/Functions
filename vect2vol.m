function volout = vect2vol(vect,atlas)

volout = atlas;
for i = 1:length(vect)
   volout(find(atlas==i)) = vect(i); 
end