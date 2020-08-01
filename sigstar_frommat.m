function sigstar_frommat(xlocs,pmat,fsize)

if nargin < 3
   fsize = 18; 
end

comp = {}; pplot = [];
for i = 1:size(pmat,1)
    for ii = 1:(i-1)
        if pmat(i,ii) < 0.05
            comp = [comp {[xlocs(i) xlocs(ii)]}];
            pplot = [pplot pmat(i,ii)];
        end
    end
end

sigstar(comp,pplot,0,fsize)