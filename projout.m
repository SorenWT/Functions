function xdef = projout(x,dir)

x = vert(x); dir = vert(dir);

xdef = x-dir*(x'*dir./(norm(dir)^2));