function [pxx,f]=psdplot(data,fs,range)
% Quick function for plotting power spectra

if size(data,2) == 1
    data = data';
end

if nargin > 2
    for c = 1:size(data,1)
        [pxx(c,:),f(c,:)] = pwelch(data(c,:),[],[],[],fs); %want 3 cycles of lowest frequency in window
        frange = intersect(find(f(c,:) > range(1)),find(f(c,:) < range(2)));
        loglog(f(c,frange),pxx(c,frange),'LineWidth',3)
        hold on
    end
    xlabel('Log Frequency (Hz)')
    ylabel('Log PSD')
else
    for c = 1:size(data,1)
        [pxx(c,:),f(c,:)] = pwelch(data(c,:),[],[],[],fs);
        loglog(f(c,2:end),pxx(c,2:end),'LineWidth',3)
        hold on
    end
    xlabel('Log Frequency (Hz)')
    ylabel('Log PSD')
end




