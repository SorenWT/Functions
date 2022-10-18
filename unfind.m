function lgc = unfind(idx, N)
    %Go from indicies into logical (for vectors only)
    lgc = false(N,1);
    lgc(idx) = true;
end