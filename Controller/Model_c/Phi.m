function [ phi ] = Phi( nu )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
v=nu(2);
IvI=abs(v);
r=nu(3);
IrI=abs(r);

phi=[0 0 0 0 0 0 0 0;
     -IrI*v -r -IvI*r -IrI*r 0 0 0 0;
     0 0 0 0 -IrI*v -r -IvI*r -IrI*r];

end

