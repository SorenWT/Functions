def fooof_timefreq(in_file,out_file):
    import fooof
    import numpy as np
    from scipy.io import loadmat,savemat
    
    # Load the file - should be outputs from a Fieldtrip time-frequency analysis
    data = loadmat(in_file)
    sets = data['fooof_settings']
    data = data['freqdata']
    
    freqs = np.squeeze(data['freq'])
    freqs = freqs[()]
    freqs = np.squeeze(freqs)
    psd = np.squeeze(data['powspctrm'])
    psd = psd[()]
    data=[] #save memory
    
    # Initialize fooof
    #fm = fooof.FOOOF(sets['peak_width_limits'],sets['max_n_peaks'],sets['min_peak_amplitude'],sets['peak_threshold'],sets['background_mode'],sets['verbose']);
    fm = fooof.FOOOF(aperiodic_mode='knee')
    
    siz = psd.shape
    
    allmdls = np.empty(shape=[siz[0],siz[1],siz[3]],dtype=dict)
    
    for i in range(siz[0]):
    	for ii in range(siz[1]):
    		for iii in range(siz[3]): # assume psd is in dimension 3, as it usually is
    			thispsd=np.squeeze(psd[i,ii,:,iii])
    			allmdls[i,ii,iii]
    			if ~np.isnan(thispsd).any():
    				fm.fit(freqs,thispsd,[freqs[0],freqs[-1]])
    				allmdls[i,ii,iii]['aperiodic_params'] = fm.get_params('aperiodic_params')
    				allmdls[i,ii,iii]['peak_params'] = fm.get_params('peak_params')
    				allmdls[i,ii,iii]['gaussian_params'] = fm.get_params('gaussian_params')
    				allmdls[i,ii,iii]['error'] = fm.get_params('error')
    				allmdls[i,ii,iii]['r_squared'] = fm.get_params('r_squared')
    				
    		
    
    mdldict = {'params',allmdls}
    savemat(out_file,mdldict)
    
    return
