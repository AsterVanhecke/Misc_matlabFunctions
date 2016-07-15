function [XIn] = dblTformInv(tformDbl,XBase)

xBase = XBase(:,1);
yBase = XBase(:,2);
[x1,y1] = tforminv(tformDbl{1},xBase,yBase);
[xIn,yIn] = tforminv(tformDbl{2},x1,y1);
XIn =[xIn,yIn];

