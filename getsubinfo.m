function [info] = getsubinfo(fiffile)
% fiffile should be the full file path


outfile = extractBefore(fiffile,'.fif');
outfile = outfile{1};
outfile = [outfile '_info.json'];

pyscript = fopen([fiffile '_pyscript.py'],'w');
fprintf(pyscript,'import sys \n')
fprintf(pyscript,'sys.path.insert(0, ''/home/soren/Documents/MATLAB/Functions'') \n')
fprintf(pyscript,'sys.path.insert(0, ''/home/sorenwt/projects/def-gnorthof/sorenwt/MATLAB/Functions'') \n')
fprintf(pyscript,'from getsubinfo import getsubinfo \n')
fprintf(pyscript,['getsubinfo(''' fiffile...
    ''',''' outfile ''')'])
system(['python ' fiffile '_pyscript.py'])
system(['rm ' fiffile '_pyscript.py'])

info = jsonread(outfile);
system(['rm ' outfile])


