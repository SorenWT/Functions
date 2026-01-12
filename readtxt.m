function str = readtxt(filename)

fid = fopen(filename);
raw = fread(fid,inf);
str = char(raw');