% wildcard strcmp
%
% usage: wstrcmp(s1, s2);
%        compare s1 and s2 using wildcards
%  e.g.: wstrcmp('XY*Z', 'XYAZ');       
%   by:  justin gardner
% date:  11/11/97
function match = wstrcmp(s1, s2)

if (isempty(s1) | isempty(s2) | (length(s1) > length(s2)))
  match = 0;
  return;
end

match = 1;

for i = 1:length(s1)
  if (s1(i) ~= '*')
    if (s1(i) ~= s2(i)), match = 0; end;
  end
end
