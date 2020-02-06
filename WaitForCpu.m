function WaitForCpu(targUsage,command,varargin)
%
% 
% Inputs: 
%     targUsage (float): the target usage for the function (in percent). If usage is below this value  
%         for a certain length of time, the command will be run
%     command (string): the command to be run
%     
% Optional inputs: 
%     WaitForCpu checks a the CPU usage at regular intervals, and only runs the command if the past n 
%         checks are successful
%     historylength: the number of past successful checks needed (default = 6)
%     checkinterval: the interval (in seconds) between successive CPU checks (default = 300 seconds, 
%         ie. 5 minutes)
%     

varargin = setdefault(varargin,'historylength',6);
varargin = setdefault(varargin,'checkinterval',300);

avgusage = ones(1,EasyParse(varargin,'historylength'))*100;
checkinterval = EasyParse(varargin,'checkinterval');

while 1
    [~,usage] = system('top -bn2 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk ''{print 100 - $1}''');
    usage = str2num(usage);
    fprintf([num2str(usage(2)) newline])
    avgusage = avgusage([2:end 1]);
    avgusage(end) = usage(2);
    
    if isempty(find(avgusage > targUsage))
       eval(command);
       break;
    end
    
    pause(checkinterval) % pause before checking again
end

end
