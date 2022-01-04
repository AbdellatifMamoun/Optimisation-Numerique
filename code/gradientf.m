function [grad] = gradientf(x)
%GRADIENTF Summary of this function goes here
%   Detailed explanation goes here
x1 = 4*(x(1)+x(2)+x(3)-3)+2*(x(1)-x(2));
x2 = 4*(x(1)+x(2)+x(3)-3)-2*(x(1)-x(2))+2*(x(2)-x(3));
x3 = 4*(x(1)+x(2)+x(3)-3)-2*(x(2)-x(3));
grad = [x1;x2;x3];
end

