function [y] = gradient2(x)
%HESIC Summary of this function goes here
%   Detailed explanation goes here
y1 = -400*x(1)*(x(2)-x(1)^2)-2*(1-x(1));
y2 = 200*(x(2)-x(1)^2);
y=[y1;y2];
end
