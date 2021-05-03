function p = lmplot(X,Y,coeffs,coeffnames,yl,varargin)

argsin = varargin;

argsin = setdefault(argsin,'coeffplot','bar');

if CheckInput(argsin,'panel')
    p = EasyParse(argsin,'panel');
    pindx = EasyParse(argsin,'panelindx');
    dopanel = 1;
else
    dopanel = 0;
end

cplottype = EasyParse(argsin,'coeffplot');

if dopanel
    p(pindx{:}).pack('h',{1/2 1/2})
    
    p(pindx{:},1).select()
else
    subplot(1,2,1)
end

pred = X*coeffs;
nicecorrplot(pred,Y,{['Predicted ' yl],yl},'Plot','off')
if dopanel
    p(pindx{:},2).select()
else
    subplot(1,2,2)
end

if strcmpi(cplottype,'bar')
    bar(coeffs)
    xlabel(coeffnames)
    ylabel('Coefficient')
    FixAxes(gca,14)
elseif strcmpi(cplottype,'wordcloud')
    wordcloud(coeffnames,NormOntoRange(coeffs,[0 1]))
end
