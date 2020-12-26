function pbody = embody_json2map(res)

for i = 1:length(res)
    xd = res{i}.arrXD; yd = res{i}.arrYD; 
    pbox1l = round(res{i}.pbox1l); pbox2l = round(res{i}.pbox2l); pboxt = round(res{i}.pbox1t); % should be the same top
    
    pbox1 = zeros(524,175); pbox2 = pbox1;
    
    for q = 1:length(xd)
       if xd(q) < pbox2l 
           ycoord = round(yd(q)-pboxt); xcoord = round(xd(q)-pbox1l);
           if ycoord < 1 || ycoord > size(pbox1,1) || xcoord < 1 || xcoord > size(pbox1,2)
               warning(['Data point omitted - coordinates y=' num2str(ycoord) ' x=' num2str(xcoord)])
           else
               pbox1(ycoord,xcoord) = 1;
           end
       else
           ycoord = round(yd(q)-pboxt); xcoord = round(xd(q)-pbox2l);
           if ycoord < 1 || ycoord > size(pbox2,1) || xcoord < 1 || xcoord > size(pbox2,2)
               warning(['Data point omitted - coordinates y=' num2str(ycoord) ' x=' num2str(xcoord)])
           else
               pbox2(ycoord,xcoord) = 1;
           end
       end
    end
    
    h=fspecial('gaussian',[15 15],5); % should maybe replace this with a disk?
    pbox1=imfilter(pbox1,h); pbox2 = imfilter(pbox2,h);
    
    pbody(:,:,i) = pbox1-pbox2;
end