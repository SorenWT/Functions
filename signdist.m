function sgdist = signdist(d1,d2)

if (d1==0 && d2==0)
    sgdist = 0;
elseif (d1~=0 && d2 ~=0)
    sgdist = 0;
else
    sgdist = 1;
end