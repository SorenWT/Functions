function [pfOut] = PeakFreq_EEG_wrapper(EEG,frange)

if nargin < 2
    frange = [8 13];
end

try 
    test = fast_median(rand(10,1));
    medfunc = @fast_median;
catch
    medfunc = @median;
end

disp(' ')
disp('Computing peak frequency...')

% Filter according to Cohen (2014)
trans_width    = .15;
idealresponse  = [ 0 0 1 1 0 0 ];
filtfreqbounds = [ 0 (1-trans_width)*frange(1) frange(1) frange(2) frange(2)*(1+trans_width) EEG.srate/2 ]/(EEG.srate/2);
filt_order     = round(3*(EEG.srate/frange(1)));
filterweights  = firls(filt_order,filtfreqbounds,idealresponse);

filterdata = zeros(size(EEG.data));
for chani=1:EEG.nbchan
    filterdata(chani,:,:) = reshape( filtfilt(filterweights,1,double(reshape(EEG.data(chani,:,:),1,EEG.pnts*EEG.trials))) ,EEG.pnts,EEG.trials);
end
EEG.data = filterdata;

% use fieldtrip's hilbert transform
data = eeglab2fieldtrip(EEG,'preprocessing','none');
cfg = []; cfg.hilbert = 'complex';
filtdata = ft_preprocessing(cfg,data);

n_order = 10;
orders = linspace(10,400,n_order)/2; 
orders = round( orders/(1000/EEG.srate) );

times2saveidx = (orders(end)+1):(EEG.pnts-orders(end)); % take all data points where you have enough data for the filter

angles = unwrap(angle(filtdata.trial{1}),[],2);
fslide = diff(filtdata.fsample*angles,1,2)/(2*pi);

for oi=1:n_order
    for ti=1:length(times2saveidx)
        fslide_filt(:,oi,ti) = medfunc(fslide(:,max(times2saveidx(ti)-orders(oi),1):min(times2saveidx(ti)+orders(oi),EEG.pnts-1))');
    end
end


pfOut = mean(median(fslide_filt,2),3);
