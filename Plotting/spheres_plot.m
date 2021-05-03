function s = spheres_plot(X,Y,Z,siz,cols)

nfaces = 60;

if any(size(cols)==1)
    cols = horz(cols);
    cols = repmat(cols,length(X),1);
end

[xi,yi,zi] = sphere(nfaces);

stretchfact(1) = (max(X)-min(X)); 
stretchfact(2) = (max(Y)-min(Y));
stretchfact(3) = (max(Z)-min(Z));

for i = 1:length(X)
    xplot{i} = xi.*siz(i)*stretchfact(1)+X(i);
    yplot{i} = yi.*siz(i)*stretchfact(2)+Y(i);
    zplot{i} = zi.*siz(i)*stretchfact(3)+Z(i);
end

for i = 1:length(xplot)
    colgrad = customcolormap([0 0.5 1],[palecol(cols(i,:),0.5); cols(i,:); cols(i,:).*0.5],nfaces+1);
    colmat = repmat(permute(colgrad,[1 3 2]),1,nfaces+1);
    s(i) = surf(xplot{i},yplot{i},zplot{i},colmat,'EdgeColor','none');
    material(s(i),[0.45 0.55 0.3 20 1])
    hold on
end

lighting gouraud
axis square