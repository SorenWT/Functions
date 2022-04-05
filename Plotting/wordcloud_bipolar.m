function wc = wordcloud_bipolar(words,sizes,cmap,nhighlights)

if nargin < 3
   load('lkcmap2'); cmap = lkcmap2; 
end

if nargin < 4
   nhighlights = 8; 
end



cdata = zeros(length(sizes),3);
cdata(sizes<0,:) = repmat(cmap(round(size(cmap,1)/3),:),sum(sizes<0),1); 
cdata(sizes>0,:) = repmat(cmap(round(2*size(cmap,1)/3),:),sum(sizes>0),1);

[~,sortindx] = sort(sizes);

cdata(sortindx(1:nhighlights/2),:) = repmat(cmap(1,:),nhighlights/2,1);
cdata(sortindx((end-nhighlights/2+1):end),:) = repmat(cmap(end,:),nhighlights/2,1);


wc = wordcloud(words,abs(sizes),'Color',cdata);