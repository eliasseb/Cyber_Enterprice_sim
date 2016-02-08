function [Eta_d,Eta_dtheta,Eta_dtheta2,v_s,v_stheta,v_st,theta_r]=Guidance_elip(theta, u_d, u_ddot,h)
syms th  

Eta_dtheta=zeros(3,1);
Eta_dtheta2=zeros(3,1);



Eta_d = [5+4.5*cos((pi/180)*theta) -0.75-2.25*sin((pi/180)*theta) -pi*3/2-atan2(-4.5*sin((pi/180)*theta),-2.25*cos((pi/180)*theta))]';

Y=Eta_d_elipse(th);
dYdth=diff(Y);
dYddth=diff(Y,2);
Eta_dtheta=vpa(subs(dYdth,th,theta));
Eta_dtheta2=vpa(subs(dYddth,th,theta));
v_s=vs(u_d,theta );
YY=vs(u_d,th);
dYYdth=diff(YY);
v_stheta=vpa(subs(dYYdth,th,theta));
v_st=vs(u_ddot,theta);

theta_dot=v_s;
theta=euler2(theta_dot,theta,h); %% Integrating with euler
theta_r=theta;
end


function [ Eta_d ] = Eta_d_elipse(theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Eta_d=[5+4.5*cos((pi/180)*theta) -0.75-2.25*sin((pi/180)*theta) 3/2*pi*atan2(-4.5*sin((pi/180)*theta),-2.25*cos((pi/180)*theta))]';

end


function [ v ] = vs(u_d,theta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    v=u_d/(sqrt((4.5*pi/180*sin((pi/180)*theta)).^2+2.25*pi/180*cos((pi/180)*theta)).^2);
end


