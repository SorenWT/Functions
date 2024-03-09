function [pos,pnts,whichview] = openpose_readposture_stand(filein,extended,manualview)
% extended is either 1 or 0 (for now - may add higher numbers) defining
% whether to extract only the Posturescreen features or other ones (like
% head tilt, etc)

if nargin < 2
    extended = 0;
end

if nargin < 3
    manualview = 0;
end

datin = jsonread(filein);

for i = 1:length(datin.people)
    personconfs(i) = mean(datin.people(i).pose_keypoints_2d([3:3:length(datin.people(i).pose_keypoints_2d)]));
end

if length(datin.people) ==0
   warning('No people found in image')
   pos = NaN; pnts = NaN; whichview = NaN; return
end


if length(datin.people) > 1 && (max(personconfs) < 0.5 || sum(personconfs>0.3) > 1)
    %warning(['More than one person found in ' filein ': check to make sure the right person was selected'])
    for i = 1:length(datin.people)
        subplot(1,length(datin.people),i);
        tmpx = datin.people(i).pose_keypoints_2d([1:3:length(datin.people(i).pose_keypoints_2d)]);
        tmpy = datin.people(i).pose_keypoints_2d([2:3:length(datin.people(i).pose_keypoints_2d)]);
        tmpx(tmpx==0) = NaN; tmpy(tmpy==0) = NaN;
        scatter(tmpx,tmpy);
        set(gca,'YDir','reverse')
        set(gcf,'units','normalized')
        figpos = get(gcf,'position')
        set(gcf,'position',[0 1-figpos(4) figpos(3) figpos(4)])
        Normalize_Ylim(gcf,0)
        axis equal
        title(['Person ' num2str(i) ': mean confidence = ' num2str(mean(datin.people(i).pose_keypoints_2d([3:3:length(datin.people(i).pose_keypoints_2d)])))])
    end
    whichperson = input(['File ' filein ': ' newline 'which person is the relevant one? Please input a number from 1 to ' num2str(length(datin.people)) '.'],'s');
    whichperson = str2num(whichperson);
    if isnumeric(whichperson)
        datin.people = datin.people(whichperson);
    end
    close(gcf)
    overwrite = input(['Selected person ' num2str(whichperson) ': overwrite JSON file? (y/n)'],'s');
    if strcmpi(overwrite,'y')
        textout = jsonencode(datin);
        writefile(filein,textout);
        disp(['File ' filein ' overwritten: irrelevant person removed']);
    end
else
    [~,whichperson] = max(personconfs);
    datin.people = datin.people(whichperson);
    disp(['Person ' num2str(whichperson) ' selected: mean confidence ' num2str(personconfs(whichperson))])
end
datin = datin.people(1);
datin.pose_keypoints_2d(datin.pose_keypoints_2d==0) = NaN;

pnts = struct;
pnts.x = datin.pose_keypoints_2d([1:3:length(datin.pose_keypoints_2d)]);
pnts.y = datin.pose_keypoints_2d([2:3:length(datin.pose_keypoints_2d)]);
pnts.c = datin.pose_keypoints_2d([3:3:length(datin.pose_keypoints_2d)]);

pnts.rawpos = [vert(pnts.x) vert(pnts.y)];
pnts.pos = pnts.rawpos-min(pnts.rawpos,[],1);
normfact = 72./(max(pnts.rawpos(:,2))-min(pnts.rawpos(:,2)));
pnts.pos = pnts.pos.*normfact;

if isfield(datin,'augpts')
    datin.augpts = (datin.augpts-min(pnts.rawpos,[],1)).*normfact;
end

if isfield(datin,'whichview')
    whichview = datin.whichview;
else
    hipdist = norm(pnts.pos(10,:)-pnts.pos(13,:));
    
    % figure out what view the photo is from
    
    if any(isnan(pnts.pos(18,:))) || any(isnan(pnts.pos(19,:))) %(hipdist < 7 || isnan(hipdist)
        
        %headangle = pnts.pos(1,:)-pnts.pos(2,:);
        %headangle = rad2deg(asin(headangle(1)/norm(headangle)));
        
        if any(isnan(pnts.pos(18,:)))
            whichview = 'left';
        else
            whichview = 'right';
        end
        
        %if headangle > 0
        %    whichview = 'right';
        %elseif headangle < 0
        %    whichview = 'left';
        %end
    elseif isnan(pnts.x(1))
        whichview = 'posterior';
    elseif ~isnan(pnts.x(16)) && ~isnan(pnts.x(17))
        whichview = 'anterior';
    else
        whichview = 'unknown';
    end
    
    % figure out if standing or seated
    % use relative y distance b/w hips to knees and knees to ankles - should be small for seated
    
    avg_hipknee_y = nanmean([pnts.y(11)-pnts.y(10) pnts.y(14)-pnts.y(13)]);
    avg_kneeankle_y = nanmean([pnts.y(12)-pnts.y(11) pnts.y(15)-pnts.y(14)]);
    thigh_rat = avg_hipknee_y./avg_kneeankle_y;
    
    hip_x = abs([pnts.x(10)-pnts.x(13)]); knee_x = abs([pnts.x(11)-pnts.x(14)]);
    hipknee_rat = knee_x./hip_x;
    
    %if (thigh_rat < 0.6 || isnan(thigh_rat)) || (hipknee_rat > 1.33 && ~contains(whichview,{'left','right'}))
    %    sitstand = 'seated';
    %else
    %    sitstand = 'standing';
    %end
    
    %whichview = [whichview '_' sitstand];
    
    [path,name] = fileparts(filein);
    
    path = tokenize(path,'/');
    path{end} = 'openpose_output_img';
    path = fullfile(path{:});
    
    imgfile = fullfile(['/' path],replace(name,'keypoints','rendered.jpg'));
    if manualview
        if exist(imgfile,'file')
            I = imread(imgfile);
            imshow(I)
            set(gcf,'units','normalized')
            figpos = get(gcf,'position')
            set(gcf,'position',[0 1-figpos(4) figpos(3) figpos(4)])
            title(['Automatically determined view: ' whichview],'FontSize',16)
            
            disp(['File: ' filein])
            disp(['Determined view: ' whichview])
            acceptview = input(['Is this correct? Press enter if yes, otherwise enter the correct view: '],'s');
            
            if ~isempty(acceptview)
                whichview = acceptview;
                %             switch acceptview
                %                 case 'a'
                %                     whichview = 'anterior';
                %                 case 'p'
                %                     whichview = 'posterior';
                %                 case 'l'
                %                     whichview = 'left';
                %                 case 'r'
                %                     whichview = 'right';
                %             end
            end
            
            datin.whichview = whichview;
            
            writedata = struct; writedata.people(1) = datin;
            
            
            textout = jsonencode(writedata);
            writefile(filein,textout);
            disp(['File ' filein ' updated: view ' whichview ' noted']);
            close(gcf)
        else
            warning(['No openpose image file found at ' imgfile '. Proceeding with automatically determined view "' whichview '".'])
        end
    end
    
end

% version 1: get all the stuff that Posturescreen gets to validate
pos = opose_calcangles(pnts,whichview,extended,datin);

pnts.whichview = whichview; % save for posterity