function afni_nii2img(niftiin)

if niftiin(1) ~= filesep
   niftiin = fullfile(pwd,niftiin); 
end
ext = extractAfter(niftiin,'.');

command = ['3dcalc -a ' niftiin ' -expr ''a'' -prefix ' replace(niftiin,ext,'img')];

setenv('PATH','/usr/local/fsl/bin:/anaconda3/bin:/Library/Frameworks/Python.framework/Versions/3.6/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
path = getenv('PATH');
setenv('PATH',[path ':/Users/Soren/abinexport PATH=$PATH:/Users/Soren/abin']);
system(['afni' newline newline command])
