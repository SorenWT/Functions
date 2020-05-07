def getsubinfo(fif_file,out_file):
    import mne
    import json
    
    # Load the file
    info = mne.io.read_info(raw_file)
    #dat = mne.io.read_raw_fif(fif_file)
    
    subinfo = info['subject_info'] 
    measdate = info['meas_date']
    subinfo['measdate'] = measdate.__str__()
    

    # Write to disk
    with open(out_file,'w') as f:
        json.dump(subinfo,f)
    
    return