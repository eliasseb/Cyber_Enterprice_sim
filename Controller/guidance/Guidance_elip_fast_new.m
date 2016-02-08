function [Eta_t,Eta_ttheta,Eta_ttheta2,v_s,v_stheta,v_st,theta_r]=Guidance_elip_fast_new(Eta_t, u_t, u_tdot,h)

persistent k theta
if (isempty(k))
    k=0;
    theta=0;
    
end
 angle_tan=-pi*3/2-atan2(-4.5*sin(theta),-4.5*cos(theta));   
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
 
 Eta_t = [5+4.5*cos(theta) -4.5*sin(theta) angle_new]';

Eta_ttheta=[-4.5*sin(theta) -4.5*cos(theta)  -1]';
Eta_ttheta2=-[-4.5*cos(theta) 4.5*sin(theta)  0]';
v_s=u_t/4.5;
v_stheta=0;
v_st=u_tdot/4.5;
theta_dot=v_s;
theta=euler2(theta_dot,theta,h); %% Integrating with euler
theta_r=theta;
end





