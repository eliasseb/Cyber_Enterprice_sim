function [ Eta_d ] = Eta_d_elipse(theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


Eta_d=[5+4.5*cos((pi/180)*theta) -0.75-2.25*sin((pi/180)*theta) -pi*3/2-atan2(-4.5*sin((pi/180)*theta),-2.25*cos((pi/180)*theta))];

end

