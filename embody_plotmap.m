function embody_plotmap(mapmat,mask,base)

if ~exist('mask','var') || isempty(mask)
    %mask = abs(mapmat)>prctile(reshape(abs(mapmat(abs(mapmat)>0)),[],1),1);
    mask = imread('~/Desktop/masters/embody-test/embody/matlab/mask.png');
    mask = [zeros(522,2) mask zeros(522,2)];
    mask = [zeros(1,175); mask; zeros(1,175)];
    %inmask = find(mask > 128);
end

if ~exist('base','var') || isempty(base)
   base = uint8(imread('base.png')); 
end

base2=base(10:531,33:203,:); % single image base


load('hotcoldmap');

M=max(abs(mapmat(:))); % max range for colorbar


imagesc(base2);
axis('off');
set(gcf,'Color',[1 1 1]);
hold on;
over2=mapmat;
fh=imagesc(over2,[-M,M]);
axis('off');
axis equal
colormap(hotcoldmap);
set(fh,'AlphaData',mask)