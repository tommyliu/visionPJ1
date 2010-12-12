function c=drawcircle(y,x,size,colour)

hold on;

t = 0:pi/50:2*pi;    
plot(size/2*sin(t)+x,size/2*cos(t)+y,colour);

hold off;
