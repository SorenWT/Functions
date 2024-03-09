function wc = wordcloud_bipolar(words,sizes,cmap,varargin)

if nargin < 3
    load('lkcmap2'); cmap = lkcmap2;
end

argsin = varargin;
if isempty(argsin)
    argsin = setdefault(argsin,'top',8);
    argsin = setdefault(argsin,'highlight','on');
end


cdata = zeros(length(sizes),3);
cdata(sizes<0,:) = repmat(cmap(round(size(cmap,1)/3),:),sum(sizes<0),1);
cdata(sizes>0,:) = repmat(cmap(round(2*size(cmap,1)/3),:),sum(sizes>0),1);


if EasyParse(argsin,'highlight','on')
    if CheckInput(argsin,'top')
        nhighlights = EasyParse(argsin,'top');
        [~,sortindx] = sort(sizes);
        
        cdata(sortindx(1:nhighlights/2),:) = repmat(cmap(1,:),nhighlights/2,1);
        cdata(sortindx((end-nhighlights/2+1):end),:) = repmat(cmap(end,:),nhighlights/2,1);
    elseif CheckInput(argsin,'specific')
        thesehighlights = EasyParse(argsin,'specific');
        if ~islogical(thesehighlights)
            thesehighlights = unfind(thesehighlights,size(cdata,1));
        end
        
        cdata(logical(thesehighlights.*(sizes<0)),:) = repmat(cmap(1,:),sum(thesehighlights.*(sizes<0)),1);
        cdata(logical(thesehighlights.*(sizes>0)),:) = repmat(cmap(end,:),sum(thesehighlights.*(sizes>0)),1);
    end
else
    cdata(sizes<0,:) = repmat(cmap(1,:),sum(sizes<0),1);
    cdata(sizes>0,:) = repmat(cmap(end,:),sum(sizes>0),1);
end


wc = wordcloud(words,abs(sizes),'Color',cdata);