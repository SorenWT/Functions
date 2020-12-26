function img = base64saveimg(fname,stringin)

%base64 = org.apache.commons.codec.binary.Base64;
%imgdat = base64.decode(uint8(stringin));

imgdat = base64decode(stringin);

f = fopen(fname,'wb');
fwrite(f,imgdat,'uint8');
fclose(f);

if nargout > 0
   img = imread(fname); 
end

