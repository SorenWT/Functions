function f = plotregline(a,b,varargin)

argsin = varargin;
argsin = setdefault(argsin,'plotCI','off');
argsin = setdefault(argsin,'color',[1 0 0]);
argsin = setdefault(argsin,'xrange',[min(a) max(a)]);

xrange = EasyParse(argsin,'xrange');
argsin = removeargs(argsin,{'xrange'});

mdl = fitlm(a,b);
B = mdl.Coefficients.Estimate;

if EasyParse(argsin,'plotCI','on')
    B_ci = coefCI(mdl); 
end

%B = regress(b,[ones(size(a)) a]);

if CheckInput(argsin,'plotCI') && EasyParse(argsin,'plotCI','on')
    hold on
    
    
    B_inci(1,:) = repmat(B(1),1,100);
    %B_inci(1,:) = linspace(B_ci(1,1),B_ci(1,2),100);
    B_inci(2,:) = linspace(B_ci(2,1),B_ci(2,2),100);
    
    for i = 1:100
       reglines(:,i) = B_inci(1,i)+B_inci(2,i)*linspace(xrange(1),xrange(2),1000);
    end
    
    curve_upper = max(reglines,[],2); curve_lower = min(reglines,[],2); 
    
    clr = EasyParse(argsin,'color');
    patch([linspace(xrange(1),xrange(2),1000) fliplr(linspace(xrange(1),xrange(2),1000))],...
        [curve_upper' fliplr(curve_lower')],clr,'EdgeColor','none','FaceAlpha',0.1);
end

argsin = removeargs(argsin,{'plotCI'});
f = plot(linspace(xrange(1),xrange(2),1000),B(1)+B(2)*linspace(xrange(1),xrange(2),1000),argsin{:});
