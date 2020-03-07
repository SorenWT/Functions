function colorGradient = createColorGradient(startColor,endColor,nColors)
% createColorGradient
% Gradient between two colors 
%
% colorGradient = createColorGradient(startColor,endColor,nColors) returns
% an RGB color gradient of length nColors from startColor to endColor.
% startColor and endColor should be RGB values with intensities in the
% range [0 1].

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

colorGradient = zeros(nColors,3);
for i = 1 : 3
    if endColor(i)~=startColor(i)
        colorGradient(:,i) = startColor(i) : (endColor(i)-startColor(i))/(nColors-1) : endColor(i);
    elseif endColor(i)==startColor(i)
        colorGradient(:,i) = startColor(i);
    end
end