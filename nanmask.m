function [mask] = nanmask(logicin)

mask = double(logicin); mask(mask==0) = NaN;