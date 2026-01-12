function [prcout] = prctile_swt(X,prc)

X = X(~isnan(X)); 
prcout = prctile(X,prc);