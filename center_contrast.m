function controut = center_contrast(contrin);

contrin = horz(contrin);

m = nanmean(contrin,2);
contrin(contrin~=0) = contrin(contrin~=0)-m;
controut = contrin;