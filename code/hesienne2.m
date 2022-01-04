function [y] = hesienne2(x)
%HESIENNE2 Summary of this function goes here
%   Detailed explanation goes here
y11 = -400*(x(2)-3*x(1)^2)+2;
y12=-400*x(1);
y21=-400*x(1);
y22=200;
y=[y11 y12;y21 y22];
end

