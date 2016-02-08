function [ Eta_d ] = Eta_d_line(theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
eta_d_0 = [ 9.5000,  -0.7500, pi]';


Eta_d=eta_d_0*(1-theta)+ theta*[0 0 pi-atan2(0.75,9.5)]';

end

