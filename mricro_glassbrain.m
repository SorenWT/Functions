function mricro_glassbrain(dat,hdr)

niftiwrite(dat,fullfile(pwd,'tmp.nii'),hdr);

file2load = fullfile(pwd,'tmp.nii');

% command = ['/Applications/MRIcroGL.app/Contents/MacOS/MRIcroGL ' ...
%     '-std -dr 2000 6000 ' char(39) file2load char(39) ' -dr ' num2str(max(max(max(dat)))) ' ' num2str(min(min(min(dat(find(dat~=0))))))];
% 
% system(command)

pyscript = fopen('mricro_pyscript.py','w');
fprintf(pyscript,'import gl \n')
fprintf(pyscript,'gl.resetdefaults() \n')
fprintf(pyscript,'gl.loadimage(''spm152'') \n')
fprintf(pyscript,['gl.overlayload(' char(39) file2load char(39) ') \n'])
fprintf(pyscript,'gl.shadername(''Glass'') \n')
fprintf(pyscript,'gl.backcolor(255,255,255)')


system(['/Applications/MRIcroGL.app/Contents/MacOS/MRIcroGL mricro_pyscript.py'])
system('rm mricro_pyscript.py')
system('rm tmp.nii')