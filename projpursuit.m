function [w,projcomps,wm,dwm,resid] = projpursuit(data,minfunc,varargin)
% data should be in the format variables x observations (like EEG data -
% i.e. if observations are subjects, they should be rows, not columns!)
% minfunc is the function to be minimized

argsin = varargin;

argsin = setdefault(argsin,'functhresh',NaN);
argsin = setdefault(argsin,'itermax',10);
argsin = setdefault(argsin,'do_whiten',0);
argsin = setdefault(argsin,'UseParallel',false);

itermax = EasyParse(argsin,'itermax');
functhresh = EasyParse(argsin,'functhresh');
do_whiten = EasyParse(argsin,'do_whiten');
do_parallel = EasyParse(argsin,'UseParallel');

% prewhiten
%data = data-mean(data,1);
if do_whiten
    [E, D] = pcamat(data);
    [wdata, wm, dwm] = whitenv(data, E, D);
    wdataorig = wdata;
else
    wm = eye(size(data,1)); dwm = eye(size(data,1));
    wdata = data;
    wdataorig = wdata;
end

nvars = size(wdata,1);


tic
w(1,:) = ga(@(west)minfunc(wdata,west),nvars,[],[],[],[],-ones(nvars,1),ones(nvars,1),@normcond,struct('UseVectorized',true,'MaxTime',300,'UseParallel',do_parallel));
w(1,:) = w(1,:)./norm(w(1,:)); % shouldn't need to do this now but just in case
toc

i = 1;

gof(i) = minfunc(wdata,w(i,:));

%thiscomp = ker(w(i,:),wdata);
thiscomp = w(i,:)*wdata;

wdata = wdata-w(i,:)'*(w(i,:)*wdata);

while 1
    i = i+1;
    tic
    w(i,:) = ga(@(west)minfunc(wdata,west),nvars,[],[],w(1:i-1,:),zeros(i-1,1),-ones(nvars,1),ones(nvars,1),@normcond,struct('UseVectorized',true,'MaxTime',300));
    if norm(w) ~= 0
        w(i,:) = w(i,:)./norm(w(i,:));
    else
        error('Weight vector converged towards zero.')
    end
    toc
    
    [gof(i)] = minfunc(wdata,w(i,:));
    
    thiscomp = w(i,:)*wdata;
    
    if(1-gof(i)) < functhresh
        w(i,:) = []; gof(i) = [];
        break;
    end
    
    if i == itermax
        break;
    end
    
    wdata = wdata-w(i,:)'*(w(i,:)*wdata);
end

gof = 1-gof;

projcomps = w*wm*data;
%fracspace = dwm*w'*projcomps;

resid = data;
for i = 1:size(w,1)
    resid = resid-dwm*w(i,:)'*w(i,:)*wm*data;
end

w = dwm*w';

%fracspace = w*origdata;
%oscispace = origdata-fracspace;

end

function [c,ceq] = normcond(w)

c = [];
ceq = norm(w)-1; % unit norm constraint

end

%
% function [y] = minfunc(x,w,frange,fs)
%
% if istable(w)
%     w = horz(w{:,:});
% end
% w = horz(w);
% %w = vert(w);
%
% %proj = ker(w,x);
%
% proj = w*x;
%
% [pxx,f] = pwelch(proj,[],[],[],fs);
%
% slope_index = (f > frange(1) & f < frange(2));
%
% linfreq = log10(f(slope_index));
% fitdata = log10(pxx(slope_index));
%
% [~,~,~,~,y] = regress(fitdata,[ones(size(linfreq)),linfreq]);
% %y = y(1);
% y = 1-y(1);
% %yorig = y(1);
% %y = 0.5-abs(0.5-yorig);
% end

