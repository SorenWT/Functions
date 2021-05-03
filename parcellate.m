function [parcvect,parcimage] = parcellate(image,atlas)
% this function assumes that zeros in the atlas are not relevant - only
% positive integers are regions

regs = unique(atlas(atlas>0));

parcvect = zeros(size(regs));
for i = 1:length(regs)
    parcvect(i) = nanmean(image(atlas==regs(i)));
    parcimage(atlas==regs(i)) = parcvect(i);
end