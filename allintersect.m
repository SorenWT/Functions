function C = allintersect(varargin)

argsin = varargin;

C = varargin{1};
for i = 2:length(varargin)
    C = intersect(C,varargin{i});
end