function [ylim] = stdshade(F,amatrix,acolor,alpha,tsdimn,method,smth)
% usage: stdshade(F,amatrix,acolor,alpha,tsdimn,method,smth)
% plot mean and sem/std coming from a matrix of data.
%sem/std is shown as shading.
% - tsdimn specifies the dimension over which to plot the data (default = 2 -
% rows are observations)
% - acolor defines the used color (default is red)
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23

%modified by SWT on 19/03/2019

%modified by SWT on 27/06/2019

if tsdimn == 1
    amatrix = amatrix';
end

if ~exist('method','var')
    method = 'sem';
end

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r';
end

if exist('F','var')==0 || isempty(F);
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1;
end

if ne(size(F,1),1)
    F=F';
end

if strcmpi(method,'mad')
    amean = smooth(nanmedian(amatrix),smth)';
else
    amean=smooth(nanmean(amatrix),smth)';
end

switch method
    case 'std'
        astd = nanstd(amatrix); % to get std shading
    case 'sem'
        astd = nanstd(amatrix)/sqrt(size(amatrix,1)-1); % to get sem shading
    case 'mad'
        astd = mad(amatrix);
    case 'logstd' % std shading when log-scaling data
        astd = nanstd(log10(amatrix));
    case 'paramci'
        astd = 1.96*nanstd(amatrix)/sqrt(size(amatrix,1)-1);
    case 'prctileci'
        astd(1,:) = prctile(amatrix,2.5,1);
        astd(2,:) = prctile(amatrix,97.5,1);
end

if exist('alpha','var')==0 || isempty(alpha)
    if strcmpi(method,'prctileci')
        fill([F fliplr(F)],[astd(2,:) fliplr(astd(1,:))],acolor,'linestyle','none','HandleVisibility','off');
        acolor='k';
    elseif ~strcmpi(method,'logstd')
        fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none','HandleVisibility','off');
        acolor='k';
    else
        fill([F fliplr(F)],[amean.*astd fliplr(amean./astd)],acolor,'linestyle','none','HandleVisibility','off');
        acolor='k';
    end
else
    if strcmpi(method,'prctileci')
        fill([F fliplr(F)],[astd(2,:) fliplr(astd(1,:))],acolor,'FaceAlpha',alpha,'linestyle','none','HandleVisibility','off');
    elseif ~strcmpi(method,'logstd')
        fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none','HandleVisibility','off');
    else
            fill([F fliplr(F)],[amean.*astd fliplr(amean./astd)],acolor, 'FaceAlpha', alpha,'linestyle','none','HandleVisibility','off');
    end
end

ylim = [min([amean+astd amean-astd]),max([amean+astd amean-astd])];

if ishold==0
    check=true; else check=false;
end

hold on;
if ischar(acolor)
    plot(F,amean,acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line
else
    plot(F,amean,'Color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line
end


if check
    hold off;
end

end



