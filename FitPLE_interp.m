function [p,ple] = FitPLE_interp(freq,pow)

tmp = log10(freq);
pow = log10(pow);
lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
pow = interp1(tmp,pow,lintmp);
p = polyfit(lintmp,pow,1);
ple = -p(1);