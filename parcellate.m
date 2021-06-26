function [parcvect,parcimage] = parcellate(image,atlas,parcelmethod)
% this function assumes that zeros in the atlas are not relevant - only
% positive integers are regions

if nargin < 3
    parcelmethod = 'mean';
end

regs = unique(atlas(atlas>0));

parcvect = zeros(size(regs));
for i = 1:length(regs)
    switch parcelmethod
        case 'mean'
            parcvect(i) = nanmean(image(atlas==regs(i)));
        case 'median'
            parcvect(i) = nanmedian(image(atlas==regs(i)));
        case 'peak'
            tmpimg = image(atlas==regs(i));
            signval = sign(tmpimg(abs(tmpimg)==max(abs(tmpimg))));
            if isempty(signval)
               signval = 1; 
            end
            signval = unique(signval);
            if length(unique(signval)) > 1
                warning('Equal magnitude peak positive and negative activations - choosing the one with most activity')
                if sign(mean(tmpimg)) > 0
                   signval = 1;
                else
                    signval = -1;
                end
            end
            if isempty(signval)
                error('test')
            end
            parcvect(i) = max(abs(tmpimg)).*signval;
    end
    parcimage(atlas==regs(i)) = parcvect(i);
end