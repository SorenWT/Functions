function xout = interactions(xin)

for i = 1:size(xin,2)
   for ii = 1:(i-1)
       xout{i,ii} = xin(:,i).*xin(:,ii);
   end
end
xout = reshape(xout,[],1);
xout = cat(2,xout{:});