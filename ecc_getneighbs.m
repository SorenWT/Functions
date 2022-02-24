function neighbs = ecc_getneighbs(datasetinfo)


switch datasetinfo.atlasname
        case 'aal'
            tissue = datasetinfo.atlas.tissue;
            tissue(find(tissue == 0)) = NaN;
            
            neighbs = struct;
            for c = 1:length(datasetinfo.atlas.tissuelabel)
                neighbs(c).label = datasetinfo.atlas.tissuelabel{c};
                neighbs(c).neighblabel = cell(1,1);
            end
            
            for c = 1:3
                edges{c} = diff(tissue,1,c);
                edges{c} = edges{c} > 0;
                boundaryindx = find(edges{c});
                for cc = 1:length(boundaryindx)
                    [i1,i2,i3] = ind2sub(size(edges{c}),boundaryindx(cc));
                    indx = {i1 i2 i3};
                    indx{c} = indx{c}+1;
                    roi1val = tissue(i1,i2,i3);
                    roi2val = tissue(indx{:});
                    neighbs(roi1val).neighblabel = [neighbs(roi1val).neighblabel neighbs(roi2val).label];
                    neighbs(roi2val).neighblabel = [neighbs(roi2val).neighblabel neighbs(roi1val).label];
                end
            end
            
            for c = 1:length(neighbs)
                neighbs(c).neighblabel(1) = [];
                neighbs(c).neighblabel = unique(neighbs(c).neighblabel);
            end
        case {'mmp','yeo'}
            
            if strcmpi(datasetinfo.atlasname,'yeo')
                datasetinfo.atlas.parcellationlabel = cellstr(num2str([1:8004]'));
                datasetinfo.atlas.parcellation = datasetinfo.atlas.parcels;
            end
            
            pos = datasetinfo.atlas.pos;
            tri = datasetinfo.atlas.tri;
            
            vox_neighbs = cell(1,length(pos));
            
            for c = 1:length(pos)
                inds = find(tri == c);
                for cc = 1:length(inds)
                    [i1,i2] = ind2sub(size(tri),inds(cc));
                    vox_neighbs{c} = [vox_neighbs{c} datasetinfo.atlas.tri(i1,except(1:3,i2))];
                end
                vox_neighbs{c} = unique(vox_neighbs{c});
            end
            
            reg_neighbs = cell(1,length(datasetinfo.atlas.parcellationlabel));
            for c = 1:length(datasetinfo.atlas.parcellationlabel)
                reg_neighbs{c} = cat(2,vox_neighbs{find(datasetinfo.atlas.parcellation == c)});
                for cc = 1:length(reg_neighbs{c})
                    reg_neighbs{c}(cc) = datasetinfo.atlas.parcellation(reg_neighbs{c}(cc));
                end
                reg_neighbs{c} = unique(reg_neighbs{c});
                reg_neighbs{c}(find(reg_neighbs{c} == c)) = [];
            end
            
            neighbs = struct;
            for c = 1:length(reg_neighbs)
                neighbs(c).label = datasetinfo.atlas.parcellationlabel{c};
                neighbs(c).neighblabel = {datasetinfo.atlas.parcellationlabel{reg_neighbs{c}}};
            end
            
            vox_neighbs = [];
            reg_neighbs = [];
            pos = [];
            tri = [];
    end