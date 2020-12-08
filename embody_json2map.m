function embody_json2map(res)

for i = 1:length(res)
    xd = res{i}.arrXD; yd = res{i}.arrYD; 
    pbox1l = round(res{i}.pbox1l); pbox2l = round(res{i}.pbox2l); pboxt = round(res{i}.pbox1t); % should be the same top
    
    pbox1 = zeros(524,175); pbox2 = pbox1;
    
    for q = 1:length(xd)
       if xd(q) < pbox2l 
           pbox1(yd(q)-pboxt,xd(q)-pbox1l) = 1;
       else
           pbox2(yd(q)-pboxt,xd(q)-pbox2l) = 1;
       end
    end
    
    h=fspecial('gaussian',[15 15],5); % should maybe replace this with a disk?
    pbox1=imfilter(pbox1,h); pbox2 = imfilter(pbox2,h);
    
    pbody(:,:,i) = pbox1-pbox2;
end