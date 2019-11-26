function cvout = CV(data,dim)

cvout = std(data,[],dim)./mean(data,dim);