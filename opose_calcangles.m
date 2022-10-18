function pos = opose_calcangles(pnts,whichview,extended,datin)

pos = table;

if contains(whichview,'_')
    whichview = extractBefore(whichview,'_');
    %whichview = whichview{1};
end

switch whichview
    case 'anterior'
        pos.Head_AnteriorAngulations = rad2deg(asin((pnts.pos(19,2)-pnts.pos(18,2))./norm(pnts.pos(19,:)-pnts.pos(18,:)))); % changed from posturescreen version (measure between eyes) to measuring between ears for consistency with posterior view
        pos.Shoulder_AnteriorAngulations = rad2deg(asin((pnts.pos(6,2)-pnts.pos(3,2))./norm(pnts.pos(6,:)-pnts.pos(3,:))));
        pos.('Hips/Pelvis_AnteriorAngulations') = rad2deg(asin((pnts.pos(13,2)-pnts.pos(10,2))./norm(pnts.pos(13,:)-pnts.pos(10,:))));
        pos.Head_AnteriorAngulationsX = rad2deg(asin((pnts.pos(1,1)-pnts.pos(2,1))./norm(pnts.pos(2,:)-pnts.pos(1,:))));
        pos.Shoulder_AnteriorAngulationsX = rad2deg(asin((pnts.pos(2,1)-pnts.pos(9,1))./norm(pnts.pos(9,:)-pnts.pos(2,:))));
        pos.Ribcage_AnteriorAngulationsX = pos.Shoulder_AnteriorAngulationsX;
        pos.('Hips/Pelvis_AnteriorAngulationsX') = rad2deg(asin((mean(pnts.pos([10 13],1),1)-mean(pnts.pos([12 14],1),1))./...
            norm(mean(pnts.pos([10 13],:),1)-mean(pnts.pos([12 14],:),1))));
        if extended
            pos.Stance_AnteriorAngulations = rad2deg((asin((pnts.pos(10,1)-pnts.pos(12,1))./norm(pnts.pos(12,:)-pnts.pos(10,:)))+...
                asin((pnts.pos(15,1)-pnts.pos(13,1))./norm(pnts.pos(15,:)-pnts.pos(13,:))))/2);
            pos.Elbows_AnteriorAngulations = rad2deg((asin((pnts.pos(3,1)-pnts.pos(4,1))./norm(pnts.pos(3,:)-pnts.pos(4,:)))+...
                asin((pnts.pos(7,1)-pnts.pos(6,1))./norm(pnts.pos(7,:)-pnts.pos(6,:))))/2);
            %pos.Feet_AnteriorAngulations = rad2deg((asin((pnts.pos(12,1)-pnts.pos(23,1))./norm(pnts.pos(12,:)-pnts.pos(23,:)))+...
            %    asin((pnts.pos(20,1)-pnts.pos(15,1))./norm(pnts.pos(20,:)-pnts.pos(15,:))))/2);
        end
    case 'posterior'
        % posterior is the same, but flip the sign of the angles
        pos.Head_PosteriorAngulations = rad2deg(asin((pnts.pos(19,2)-pnts.pos(18,2))./norm(pnts.pos(19,:)-pnts.pos(18,:))));
        pos.Shoulder_PosteriorAngulations = rad2deg(asin((pnts.pos(6,2)-pnts.pos(3,2))./norm(pnts.pos(6,:)-pnts.pos(3,:))));
        pos.('Hips/Pelvis_PosteriorAngulations') = rad2deg(asin((pnts.pos(13,2)-pnts.pos(10,2))./norm(pnts.pos(13,:)-pnts.pos(10,:))));
        pos.Head_PosteriorAngulationsX = -rad2deg(asin((mean(pnts.pos([18 19],1),1)-pnts.pos(2,1))./...
            norm(pnts.pos(2,:)-mean(pnts.pos([18 19],:),1))));
        pos.Shoulder_PosteriorAngulationsX = -rad2deg(asin((pnts.pos(2,1)-pnts.pos(9,1))./norm(pnts.pos(9,:)-pnts.pos(2,:))));
        pos.Ribcage_PosteriorAngulationsX = pos.Shoulder_PosteriorAngulationsX;
        pos.('Hips/Pelvis_PosteriorAngulationsX') = -rad2deg(asin((mean(pnts.pos([10 13],1),1)-mean(pnts.pos([12 14],1),1))./...
            norm(mean(pnts.pos([10 13],:),1)-mean(pnts.pos([12 14],:),1))));
        if extended
            pos.Stance_PosteriorAngulations = -rad2deg((asin((pnts.pos(10,1)-pnts.pos(12,1))./norm(pnts.pos(12,:)-pnts.pos(10,:)))+...
                asin((pnts.pos(15,1)-pnts.pos(13,1))./norm(pnts.pos(15,:)-pnts.pos(13,:))))/2);
            pos.Elbows_PosteriorAngulations = -rad2deg((asin((pnts.pos(3,1)-pnts.pos(4,1))./norm(pnts.pos(3,:)-pnts.pos(4,:)))+...
                asin((pnts.pos(7,1)-pnts.pos(6,1))./norm(pnts.pos(7,:)-pnts.pos(6,:))))/2);
            %pos.Feet_AnteriorAngulations = rad2deg((asin((pnts.pos(12,1)-pnts.pos(23,1))./norm(pnts.pos(12,:)-pnts.pos(23,:)))+...
            %    asin((pnts.pos(20,1)-pnts.pos(15,1))./norm(pnts.pos(20,:)-pnts.pos(15,:))))/2);
        end
    case 'right'
        pos.Knees_LateralAngulations = rad2deg(asin((pnts.pos(11,1)-pnts.pos(12,1))./norm(pnts.pos(11,:)-pnts.pos(12,:))));
        pos.('Hips/Pelvis_LateralAngulations') = rad2deg(asin((pnts.pos(10,1)-pnts.pos(11,1))./norm(pnts.pos(10,:)-pnts.pos(11,:))));
        pos.Shoulder_LateralAngulations = rad2deg(asin((pnts.pos(3,1)-pnts.pos(10,1))./norm(pnts.pos(3,:)-pnts.pos(10,:))));
        pos.Head_LateralAngulations = rad2deg(asin((pnts.pos(18,1)-pnts.pos(3,1))./norm(pnts.pos(18,:)-pnts.pos(3,:))));
        if extended
            pos.Head_LateralAngulationsX = -rad2deg(asin((pnts.pos(1,2)-pnts.pos(18,2))./norm(pnts.pos(1,:)-pnts.pos(18,:))));
            pos.Elbows_LateralAngulations = rad2deg(asin((pnts.pos(4,1)-pnts.pos(3,1))./norm(pnts.pos(4,:)-pnts.pos(3,:))));
            if isfield(datin,'augpts')
                pos.Ribcage_LateralAngulationsLower =  rad2deg(asin((datin.augpts(3,1)-pnts.pos(10,1))./norm(datin.augpts(3,:)-pnts.pos(10,:))));
                pos.Ribcage_LateralAngulationsUpper = rad2deg(asin((pnts.pos(3,1)-datin.augpts(3,1))./norm(datin.augpts(3,:)-pnts.pos(3,:))));
                pos.Ribcage_LateralAngulationsUpperCN = rad2deg(asin((datin.augpts(2,1)-datin.augpts(3,1))./norm(datin.augpts(3,:)-datin.augpts(2,:))));
                pos.Shoulder_LateralAngulationsAC =  rad2deg(asin((datin.augpts(1,1)-pnts.pos(10,1))./norm(datin.augpts(1,:)-pnts.pos(10,:))));
                pos.Head_LateralAngulationsAC = rad2deg(asin((pnts.pos(18,1)-datin.augpts(1,1))./norm(pnts.pos(18,:)-datin.augpts(1,:))));
                pos.Chest_LateralAngulations = -rad2deg(asin((datin.augpts(2,2)-datin.augpts(1,2))./norm(datin.augpts(1,:)-datin.augpts(2,:))));
                pos.Shoulder_LateralAngulationsCN =  rad2deg(asin((datin.augpts(2,1)-pnts.pos(10,1))./norm(datin.augpts(2,:)-pnts.pos(10,:))));
                pos.Head_LateralAngulationsCN = rad2deg(asin((pnts.pos(18,1)-datin.augpts(2,1))./norm(pnts.pos(18,:)-datin.augpts(2,:))));
                pos.Shoulder_RetractionAngulation = rad2deg(atan((datin.augpts(1,2)-pnts.pos(3,2))./(datin.augpts(1,1)-pnts.pos(3,1))));
                if datin.augpts(1,1) > pnts.pos(3,1)
                    pos.Shoulder_RetractionAngulation = 180+pos.Shoulder_RetractionAngulation;
                end
            else
                pos.Ribcage_LateralAngulationsLower = NaN;
                pos.Ribcage_LateralAngulationsUpper = NaN;
                pos.Shoulder_LateralAngulationsAC =  NaN;
                pos.Head_LateralAngulationsAC = NaN;
                pos.Chest_LateralAngulations = NaN;
                pos.Shoulder_LateralAngulationsCN = NaN;
                pos.Head_LateralAngulationsCN = NaN;
            end
        end
    case 'left'
        % don't flip sign to mimic the sign issues in posturescreen
        pos.Knees_post_LateralAngulations = rad2deg(asin((pnts.pos(14,1)-pnts.pos(15,1))./norm(pnts.pos(14,:)-pnts.pos(15,:))));
        pos.('Hips/Pelvis_post_LateralAngulations') = rad2deg(asin((pnts.pos(13,1)-pnts.pos(14,1))./norm(pnts.pos(13,:)-pnts.pos(14,:))));
        pos.Shoulder_post_LateralAngulations = rad2deg(asin((pnts.pos(6,1)-pnts.pos(13,1))./norm(pnts.pos(6,:)-pnts.pos(13,:))));
        pos.Head_post_LateralAngulations = rad2deg(asin((pnts.pos(19,1)-pnts.pos(6,1))./norm(pnts.pos(19,:)-pnts.pos(6,:))));
        if extended
            pos.Head_post_LateralAngulationsX = -rad2deg(asin((pnts.pos(1,2)-pnts.pos(19,2))./norm(pnts.pos(1,:)-pnts.pos(19,:))));
            pos.Elbows_post_LateralAngulations = -rad2deg(asin((pnts.pos(7,1)-pnts.pos(6,1))./norm(pnts.pos(7,:)-pnts.pos(6,:))));
            if isfield(datin,'augpts')
                pos.Ribcage_post_LateralAngulationsLower =  -rad2deg(asin((datin.augpts(3,1)-pnts.pos(13,1))./norm(datin.augpts(3,:)-pnts.pos(13,:))));
                pos.Ribcage_post_LateralAngulationsUpper = -rad2deg(asin((pnts.pos(6,1)-datin.augpts(3,1))./norm(datin.augpts(3,:)-pnts.pos(6,:))));
                pos.Ribcage_LateralAngulationsUpperCN = -rad2deg(asin((datin.augpts(2,1)-datin.augpts(3,1))./norm(datin.augpts(3,:)-datin.augpts(2,:))));
                pos.Shoulder_post_LateralAngulationsAC =  -rad2deg(asin((datin.augpts(1,1)-pnts.pos(13,1))./norm(datin.augpts(1,:)-pnts.pos(13,:))));
                pos.Head_post_LateralAngulationsAC = -rad2deg(asin((pnts.pos(19,1)-datin.augpts(1,1))./norm(pnts.pos(19,:)-datin.augpts(1,:))));
                pos.Chest_post_LateralAngulations = -rad2deg(asin((datin.augpts(2,2)-datin.augpts(1,2))./norm(datin.augpts(1,:)-datin.augpts(2,:))));
                pos.Shoulder_post_LateralAngulationsCN =  -rad2deg(asin((datin.augpts(2,1)-pnts.pos(13,1))./norm(datin.augpts(2,:)-pnts.pos(13,:))));
                pos.Head_post_LateralAngulationsCN = -rad2deg(asin((pnts.pos(19,1)-datin.augpts(2,1))./norm(pnts.pos(19,:)-datin.augpts(2,:))));
                pos.Shoulder_post_RetractionAngulation = -rad2deg(atan((datin.augpts(1,2)-pnts.pos(6,2))./(datin.augpts(1,1)-pnts.pos(6,1))));
                if datin.augpts(1,1) < pnts.pos(6,1)
                    pos.Shoulder_post_RetractionAngulation = 180+pos.Shoulder_post_RetractionAngulation;
                end
            else
                pos.Ribcage_post_LateralAngulationsLower = NaN;
                pos.Ribcage_post_LateralAngulationsUpper = NaN;
                pos.Shoulder_post_LateralAngulationsAC =  NaN;
                pos.Head_post_LateralAngulationsAC = NaN;
                pos.Chest_post_LateralAngulations = NaN;
                pos.Shoulder_post_LateralAngulationsCN = NaN;
                pos.Head_post_LateralAngulationsCN = NaN;
            end
        end
end