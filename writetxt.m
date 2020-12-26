function writetxt(filename,textin)

f = fopen(filename,'w');
fwrite(f,textin);
fclose(f)