function P = plotci_ellipse(data2d,conf,varargin)

l = lines;

argsin = varargin; 
argsin = setdefault(argsin,'color',l(1,:));
clr = EasyParse(argsin,'color'); 
argsin = removeargs(argsin,{'color'});
argsin = setdefault(argsin,'FaceAlpha',0.25);
argsin = setdefault(argsin,'EdgeColor',clr);
argsin = setdefault(argsin,'FaceColor',clr);

gmdist = fitgmdist(data2d,1);
%gmfunc = @(x,y) arrayfun(@(x0,y0) pdf(interolowdist,[x0 y0]),x,y);
covmat = gmdist.Sigma;
[v,d] = eig(covmat);
d = diag(d); vmax = v(:,find(d==max(d))); dmax = max(d);

s = chi2inv(conf,2);

ang = atan(vmax(2,1)/vmax(1,1));

if ang < 0
    ang = ang + 2*pi;
end

R = [cos(ang) sin(ang); -sin(ang) cos(ang)];

t = linspace(0,2*pi,1000);
a = sqrt(s*d(1)); b = sqrt(s*d(2));
x = a*cos(t); y = b*sin(t);

rot = R*[horz(x); horz(y)];

x = rot(1,:)'; y = rot(2,:)';

hold on
P = patch(x+gmdist.mu(1),y+gmdist.mu(2),clr,argsin{:});
