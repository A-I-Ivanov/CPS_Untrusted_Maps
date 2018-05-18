%%This fuction tests the computation of a triangle which is tangent to the 
%%unit circle on two sides.
close all


xdoubleTic = linspace(-4, -1.1,1000);
t = linspace(-pi, pi, 100);
scaling = linspace(-5,5,1000);

for i = 1:length(xdoubleTic)
    xtic = 1/xdoubleTic(i);
    ytic = sqrt(1-xtic^2);
    
    xtic^2+ytic^2
    
    vect = [xtic - xdoubleTic(i); ytic];
    dot(vect, [xtic, ytic])
    
    m = ytic/(xtic-xdoubleTic(i));
    b = -m*xdoubleTic(i);
    
    line = vect*scaling+[xdoubleTic(i); 0;];
    plot(sin(t), cos(t), 'g');
    hold on
    
    plot(line(1,:), line(2,:), 'r');
    plot(line(1,:), -line(2,:), 'r');
    plot([xtic xtic], [ytic -ytic], 'r');
    axis([-2,2,-2,2]);
    hold off
    pause(.01);
    
end