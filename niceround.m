function rounddat = niceround(dat)

if any(size(dat)>1)
    rounddat = arrayfun(@niceround,dat);
else
    
    if abs(dat)>10
        rounddat = round(dat,1,'decimal');
    elseif abs(dat)>1
        rounddat = round(dat,2,'decimal');
    else
        rounddat = round(dat,2,'significant');
    end
end

end