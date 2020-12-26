function [a,ci]=cronbach(X,nanflag)
%Syntax: a=cronbach(X)
%_____________________
%
% Calculates the Cronbach's alpha of a data set X.
%
% a is the Cronbach's alpha.
% X is the data set. (items are columns, subjects are rows)
%
%
% Reference:
% Cronbach L J (1951): Coefficient alpha and the internal structure of
% tests. Psychometrika 16:297-333
%
%
% Alexandros Leontitsis
% Department of Education
% University of Ioannina
% Ioannina
% Greece
%
% e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
%
% June 10, 2005.
%
% Modified by SWT on 21/12/2020 to add 95% confidence intervals
% Feldt, L. S., Woodruff, D. J., & Salih, F. A. (1987). Statistical inference for coefficient alpha. Applied Psychological Measurement, 11, 93-103.



if nargin<1 | isempty(X)==1
    error('You shoud provide a data set.');
else
    % X must be a 2 dimensional matrix
    if ndims(X)~=2
        error('Invalid data set.');
    end
end

if nargin < 2
   nanflag = 0; 
end

% Calculate the number of items
if nanflag
    k=size(X,2);
    
    % Calculate the variance of the items' sum
    tmp = nansum(X'); tmp(tmp==0) = [];
    VarTotal=nanvar(tmp);
    
    % Calculate the item variance
    SumVarX=nansum(nanvar(X));
    
    % Calculate the Cronbach's alpha
    a=k/(k-1)*(VarTotal-SumVarX)/VarTotal;
else
    k=size(X,2);
    
    % Calculate the variance of the items' sum
    tmp = nansum(X'); tmp(tmp==0) = [];
    VarTotal=var(sum(X'));
    
    % Calculate the item variance
    SumVarX=sum(var(X));
    
    % Calculate the Cronbach's alpha
    a=k/(k-1)*(VarTotal-SumVarX)/VarTotal;
end

% SWT edit: added confidence intervals from Feldt et al. 1987
ci_l = 1-((1-a)*finv(1-0.05/2,size(X,1)-1,(size(X,1)-1)*(size(X,2)-1)));
ci_u = 1-((1-a)*finv(0.05/2,size(X,1)-1,(size(X,1)-1)*(size(X,2)-1)));
ci = [ci_l; ci_u];