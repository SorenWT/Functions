function z = pit(X)

lam = 1/median(X);

u = 1-exp(-lam*X);

z = normcdf(u);