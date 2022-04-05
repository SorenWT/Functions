function axpnt = figpnt2axpnt(figpnt,ax)

pos = ax.Position;
xl = ax.XLim; yl = ax.YLim;

prcpnt = [(figpnt(1)-pos(1))./pos(3) (figpnt(2)-pos(2))./pos(4)];

axpnt = [xl(1)+prcpnt(1)*range(xl) yl(1)+prcpnt(2)*range(yl)];