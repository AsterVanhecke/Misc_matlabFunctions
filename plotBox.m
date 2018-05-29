function plotBox(box,varargin)
x=zeros(1,5);
y=x;
x=x+box(1); x(2:3)=x(2:3)+box(3);
y=y+box(2); y(3:4)=y(3:4)+box(4);
if isempty(varargin)
    plot(x,y)
else
    plot(x,y,varargin{1})
end