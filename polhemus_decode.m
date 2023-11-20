function [angstbl] = polhemus_decode(polfile)

% function to get relevant angles from polhemus data

poldata = fscanf(fopen(polfile),'%f');
poldata = reshape(poldata,7,[]);
poldata([1 5:end],:) = [];
poldata(3,:) = -poldata(3,:); % flip z axis

% detect clusters
clustidx = dbscan(poldata',1.5,1);
rawpol = poldata; poldata = [];

for i = 1:length(unique(clustidx))
    poldata(:,i) = nanmean(rawpol(:,find(clustidx==i)),2);
end

if size(poldata,2) ~= 5
    scatter3(poldata(1,:),poldata(2,:),poldata(3,:)); axis equal
    
    s = input(['File ' polfile '. Assuming missing point is the nasion. Press enter if this is the case, otherwise enter n if you want to exclude this subject'],'s');
    if isempty(s)
        poldata = [NaN(3,1) poldata];
    elseif strcmpi(s,'n')
        angstbl = table; angstbl.headang = NaN; return
    end
end



% Procrustes alignment to template

template = [26.5,20,20,20,20;
    -12,-12,-19,-11.5,-4;
    33.5,33.5,25,28.5,25];

if ~any(any(isnan(poldata)))
   [~,poldata,transform] = procrustes(template,poldata,'Scaling',false);
end

poldata = array2table(poldata,'VariableNames',{'nasion','inion','right_acro','c7','left_acro'});

% calculate angles

angstbl = table;
angstbl.headang_c7 = rad2deg(atan((poldata.inion(1)-poldata.c7(1))./(poldata.inion(3)-poldata.c7(3))));
angstbl.headang_acro = rad2deg(atan((poldata.inion(1)-mean([poldata.left_acro(1) poldata.right_acro(1)]))./(poldata.inion(3)-mean([poldata.left_acro(3) poldata.right_acro(3)]))));

tmp1 = rad2deg(atan((poldata.c7(2)-poldata.right_acro(2))./(poldata.right_acro(1)-poldata.c7(1))));
tmp2 = rad2deg(atan((poldata.left_acro(2)-poldata.c7(2))./(poldata.left_acro(1)-poldata.c7(1))));
if tmp1 < 0; tmp1=180+tmp1; end;
if tmp2 < 0; tmp2=180+tmp2; end;


angstbl.shoulder_close = tmp1+tmp2;

angstbl.shoulder_raise = rad2deg(mean([atan((poldata.right_acro(3)-poldata.c7(3))./(poldata.c7(2)-poldata.right_acro(2))) ...
    atan((poldata.left_acro(3)-poldata.c7(3))./(poldata.left_acro(2)-poldata.c7(2)))]));
angstbl.head_tilt = rad2deg(atan((poldata.nasion(3)-poldata.inion(3))./(poldata.nasion(1)-poldata.inion(1))));

