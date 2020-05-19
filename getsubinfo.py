def getsubinfo(fif_file,out_file):
    import mne
    import json
    
    # Load the file
    info = mne.io.read_info(fif_file)
    #dat = mne.io.read_raw_fif(fif_file)
    
    subinfo = info['subject_info'] 
    measdate = info['meas_date']
    subinfo['measdate'] = measdate.__str__()
    subinfo.pop('height',None)
    subinfo.pop('weight',None)    

    # Write to disk
    with open(out_file,'w') as f:
        json.dump(subinfo,f)
    
    return
