function numout = neo2num(strin)

switch strin
    case 'SA'
        numout = 5;
    case 'A'
        numout = 4; 
    case 'N'
        numout = 3;
    case 'D'
        numout = 2;
    case 'SD'
        numout = 1;
end

if ~exist('numout','var')
   numout = NaN; 
end