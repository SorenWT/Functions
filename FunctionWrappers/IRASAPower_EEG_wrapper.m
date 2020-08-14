function [OsciFrac] = IRASAPower_EEG_wrapper(spec,oscifrac,frange,normalize,silent)

OsciFrac = zeros(1,size(spec.osci,2));

if nargin < 5
   silent = 0; 
end

% if nargin < 5
%     method = 'mean';
% end

if nargin < 4
    normalize = 'no';
end

if nargin < 3
    frange = [0.5 50];
end

if ~silent
    disp(' ')
    if strcmpi(oscifrac,'osci')
        disp('Computing oscillatory power...')
    else
        disp('Computing fractal power...')
    end
end

freqindex = intersect(find(spec.freq > frange(1)),find(spec.freq < frange(2)));
freqs = spec.freq(freqindex);

for c = 1:size(spec.frac,2)
    if ~silent
        fprintf([num2str(c) ' ']);
    end
    %osci = sgolayfilt(spec.osci(:,c),5,1501);
    %osci = osci(freqindex);
    %frac = sgolayfilt(spec.frac(:,c),5,1501);
    pow = spec.(oscifrac)(freqindex,c);
%     
%     if strcmpi(method,'trapz')
%         OsciFrac(c) = trapz(freqs,pow);
%     elseif strcmpi(method,'mean')
%         OsciFrac(c) = nanmean(pow);
%     end
    OsciFrac(c) = simps(freqs,pow)/numel(freqs);
    %OsciFrac(c) = (norm(pow)^2)/length(pow);
    if strcmpi(normalize,'mixd')
        OsciFrac(c) = OsciFrac(c)/(simps(freqs,specs.mixd(freqindex,c))/numel(freqs));
        %OsciFrac(c) = OsciFrac(c)/((norm(spec.mixd(freqindex,c))^2)/length(freqindex));
%         if strcmpi(method,'trapz')
%             OsciFrac(c) = OsciFrac(c)/trapz(freqs,spec.mixd(freqindex,c));
%         elseif strcmpi(method,'mean')
%             OsciFrac(c) = OsciFrac(c)/nanmean(spec.mixd(freqindex,c));
%         end
    elseif iscell(normalize)
        newfreqs = intersect(find(spec.freq > normalize{2}(1)),find(spec.freq < normalize{2}(2)));
       OsciFrac(c) = OsciFrac(c)/(simps(spec.freq(newfreqs),spec.(normalize{1})(newfreqs,c))/numel(newfreqs));
        % OsciFrac(c) = OsciFrac(c)/((norm(spec.(normalize{1})(newfreqs,c))^2)/length(newfreqs));
%         if strcmpi(method,'trapz')
%             OsciFrac(c) = OsciFrac(c)/trapz(spec.freq(newfreqs),spec.(normalize{1})(newfreqs,c));
%         elseif strcmpi(method,'mean')
%             OsciFrac(c) = OsciFrac(c)/nanmean(spec.(normalize{1})(newfreqs,c));
%         end
    end
end