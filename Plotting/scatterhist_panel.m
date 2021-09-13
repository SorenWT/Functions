function p = scatterhist_panel(X,Y,p,pindx,varargin)
% assumes for now that the histogram is always in the northeast position
% add more flexibility later

argsin = varargin;
argsin = setdefault(argsin,'loc','ne');

loc = EasyParse(argsin,'loc');

argsin = removeargs(argsin,{'loc'});


fig = figure;
s = scatterhist(X,Y,argsin{:});

switch loc
    case {'northeast','ne'}
        p(pindx{:}).pack('v',{20 80})
        p(pindx{:},1).pack('h',{80 20})
        p(pindx{:},2).pack('h',{80 20})
        
        p(pindx{:},2,1).select(s(1));
        xl = get(p(pindx{:},2,1).axis,'XLim');
        yl = get(p(pindx{:},2,1).axis,'YLim');

        % do the histograms manually because the plotting function keeps
        % messing it up
        p(pindx{:},1,1).select();
        [f,xi] = ksdensity(X);
        plot(xi,f)
        ax = gca;
        set(gca,'XLim',xl,'XTick',[])
        axis off
        ax.XAxis.Visible = 'on';
        
                p(pindx{:},2,2).select();
        [f,yi] = ksdensity(Y);
        plot(f,yi)
        ax = gca;
        set(gca,'YLim',yl,'YTick',[])
        axis off
        ax.YAxis.Visible = 'on';

        
        %p(pindx{:},2,1).select(s(1));
        %p(pindx{:},2,2).select(s(3));
end


close(fig)