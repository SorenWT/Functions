function [tortotal,tor,com] = opose_com(posein,id)

if nargin < 2
    id = NaN;
end

if ~any(isnan(id))
    scalef = opose_get_scalef(id);
else
    scalef = ones(1,10);
end

try
    [pts,allpts] = postureplot(posein,0,scalef);
catch
    tortotal = NaN; tor = NaN(1,2); com = NaN(1,3); return
end


% new connpairs specifically for COM computation
connpairs = {[1 24],[4 13],[6 15],[7 7],[9 14],[14 14],[13 16],[15 18],[16 19],[18 21]};

allpts = allpts';

for i = 1:length(connpairs)
    segpts(i,:) = (allpts(connpairs{i}(1),:)+allpts(connpairs{i}(2),:))/2; % find centers of segments
end

segmass = 0.01*[8.2,2.6,2.6,18.5,13.5,14.9,5.6,5.6,2.5,2.5]; % found online

% segmass should sum to about 1, so use it to make a weighted mean
com = sum(segpts.*segmass',1)./sum(segmass); % vector for the position of the center of mass

tor(1) = norm(com(2:3))*sin(atan(com(2)./com(3))); % inverted pendulum - torque proportional to length*sin(theta)
tor(2) = norm(com([1 3]))*sin(atan(com(1)./com(3)));
tortotal = norm(tor); % sum the x and y torques
