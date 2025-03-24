function P = plotci_ellipse(data2d,conf,varargin)

data2d = rmnan(data2d,1);

l = lines;

argsin = varargin; 
argsin = setdefault(argsin,'color',l(1,:));
clr = EasyParse(argsin,'color'); 
argsin = setdefault(argsin,'plotmeans','yes');
plotmeans = EasyParse(argsin,'plotmeans');
argsin = removeargs(argsin,{'color','plotmeans'});
argsin = setdefault(argsin,'FaceAlpha',0.25);
argsin = setdefault(argsin,'EdgeColor',clr);
argsin = setdefault(argsin,'FaceColor',clr);

gmdist = fitgmdist(data2d,1);
covmat = gmdist.Sigma;
[v,d] = eig(covmat);
d = diag(d); vmax = v(:,find(d==max(d))); dmax = max(d);

if isnumeric(conf)
    gmfunc = @(x,y) arrayfun(@(x0,y0) pdf(interolowdist,[x0 y0]),x,y);
    
    s = chi2inv(conf,2);
    
    a = sqrt(s*d(1)); b = sqrt(s*d(2));
elseif strcmpi(conf,'se')
    a = nanstd(data2d(:,2))./sqrt(numel(data2d(:,2)));
    b = nanstd(data2d(:,1))./sqrt(numel(data2d(:,1)));
elseif strcmpi(conf,'95_mean')
    b = 1.96*nanstd(data2d(:,2))./sqrt(numel(data2d(:,2)));
    a = 1.96*nanstd(data2d(:,1))./sqrt(numel(data2d(:,1)));
elseif strcmpi(conf,'sd')
    a = nanstd(data2d(:,2));
    b = nanstd(data2d(:,1));
end


   %a = a./1.96; b = b./1.96; 


ang = atan(vmax(2,1)/vmax(1,1));

if ang < 0
    ang = ang + 2*pi;
end


t = linspace(0,2*pi,1000);
x = a*cos(t); y = b*sin(t);

if isnumeric(conf)
R = [cos(ang) sin(ang); -sin(ang) cos(ang)];
rot = R*[horz(x); horz(y)];

x = rot(1,:)'; y = rot(2,:)';
end

hold on
P = patch(x+gmdist.mu(1),y+gmdist.mu(2),clr,argsin{:});
if strcmpi(plotmeans,'yes')
     scatter(gmdist.mu(1),gmdist.mu(2),160,clr,'+','linewidth',4)
end
