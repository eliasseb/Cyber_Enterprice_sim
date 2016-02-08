function R_d = Rot_dot(nu)
%% takes in the nu vector 

R_d=S([0 0 nu(3)])*Rot([0 0 nu(3)]);

end

function [ s ] = S(r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
s=[0 -r(3) r(2);
   r(3) 0 -r(1);
   -r(2) r(1) 0];

end

