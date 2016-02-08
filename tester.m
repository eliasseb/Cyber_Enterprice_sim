function [ out ] = tester( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ref=1;
 h=0.01;

out=x + (x-ref)*h;

end

