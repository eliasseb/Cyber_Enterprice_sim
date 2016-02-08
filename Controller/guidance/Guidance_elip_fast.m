function [Eta_t,Eta_ttheta,Eta_ttheta2,v_s,v_stheta,v_st,theta_r]=Guidance_elip_fast(Eta_t, u_t, u_tdot,h)

persistent k theta
if (isempty(k))
    k=0;
    theta=0;
    
end
 angle_tan=-pi*3/2-atan2(-4.5*sin((pi/180)*theta),-2.25*cos((pi/180)*theta));   
 angle_new =k*2*pi+angle_tan;
 if ((Eta_t(3)-angle_new)^2>pi^2/4)
     disp(Eta_t(3))
     disp(angle_new)
     if (angle_new>Eta_t(3))
         k=k-1;
         angle_new=k*2*pi+angle_tan;
     elseif (angle_new<Eta_t(3))
        k=k+1;
        angle_new=k*2*pi+angle_tan;
     end
 end
 
 Eta_t = [5+4.5*cos((pi/180)*theta) -0.75-2.25*sin((pi/180)*theta) angle_new]';

Eta_ttheta=[-0.0785398*sin((pi*theta)/180) -0.0392699*cos((pi*theta)/180)  -0.0349066/(4*sin((pi*theta)/180)^2+cos((pi*theta)/180)^2)]';
Eta_ttheta2=-[-0.00137078*cos((pi*theta)/180) 0.000685389*sin((pi*theta)/180)  (0.00182771*sin((pi*theta)/90))/(cos((pi*theta)/180)^2+4*sin((pi*theta)/180)^2)^2]';
v_s=u_t/(sqrt((4.5*pi/180*sin((pi/180)*theta))^2+(2.25*pi/180*cos((pi/180)*theta))^2));
v_stheta=(u_t*(0.0000807455*sin((pi*theta)/180)*cos((pi*theta)/180))/(0.0061685*sin((pi*theta)/180)^2+0.00154213*cos((pi*theta)/180)^2)^(3/2));
v_st=u_tdot/(sqrt((4.5*pi/180*sin((pi/180)*theta))^2+(2.25*pi/180*cos((pi/180)*theta))^2));
theta_dot=v_s;
theta=euler2(theta_dot,theta,h); %% Integrating with euler
theta_r=theta;
end





