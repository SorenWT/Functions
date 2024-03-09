function [jsonsurveys,pers,surveytxt] = personality_tojson(pers,fileout)
% function to export personality questionnaire structures into a surveyJS
% form

if iscell(pers)
    for i = 1:length(pers)
        [jsonsurveys{i},pers{i}] = personality_tojson(pers{i});
    end
    
    surveyjson = struct;
    surveyjson.completedHtml = '<h> Thank you for completing the personality questionnaires! The experiment will continue in a few moments. </h>';
    surveyjson.title = 'Personality module';
    surveyjson.description = 'The following form asks you a number of questions about your personality. Please answer them as honestly as you can.';
    
    for i = 1:length(jsonsurveys)
        surveyjson.pages(i).name = [jsonsurveys{i}.name '_page'];
        surveyjson.pages(i).elements = jsonsurveys(i);
    end
    
    surveytxt = jsonencode(surveyjson);
    
    surveytxt = ['var surveyJSON = ' surveytxt ';'];
    writetxt(fileout,surveytxt);
else
    jsonsurveys = struct('type','matrix','name',pers.name,'isRequired',1);
    
    if isfield(pers,'titletext')
        jsonsurveys.title = pers.titletext;
    else
        jsonsurveys.title = ['Read each statement and click the appropriate circle to its right. ' ...
            'Do not spend too much time on any statement. Answer as honestly as you can.'];
    end
    for i = 1:length(pers.respcoding)
        jsonsurveys.columns(i).value = ['Column ' num2str(i)];
        jsonsurveys.columns(i).text = pers.respcoding{i};
    end
    
    if ~isfield(pers,'prettyitems')
        pers = fix_pers_text(pers);
    end
    
    % automatically add an attention check question in each form
    tmppers = pers; thisrand = randi(length(pers.prettyitems));
    attnitem = {['Please select the "' pers.respcoding{1} '" option for this question to demonstrate your attention']};
    tmppers.prettyitems = [tmppers.prettyitems(1:thisrand-1) attnitem tmppers.prettyitems(thisrand:end)];
    for i = 1:length(tmppers.prettyitems)
        if i~=thisrand
        jsonsurveys.rows(i).value = ['Row ' num2str(i-(i>thisrand))];
        else
           jsonsurveys.rows(i).value = 'Attncheck';
        end
        jsonsurveys.rows(i).text = tmppers.prettyitems{i};
    end
    
    
end
