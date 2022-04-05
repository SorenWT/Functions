function missing = find_missing(template,data)

for i = 1:length(template)
    missing(i) = template(i)~=data(i);
    if missing(i)
       data = [data(1:i-1) NaN data(i:end)]; 
    end
end

missing = find(missing);