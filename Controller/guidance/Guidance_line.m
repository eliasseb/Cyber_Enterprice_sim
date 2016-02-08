function [Eta_d,Eta_dtheta,Eta_dtheta2,v_s,v_stheta,v_st,theta_r]=Guidance_line(Eta_d, u_d, u_ddot,h)


persistent k theta;
if (isempty(k)||isempty(theta))
    k=0;
    theta=0;
end
 eta_d_0 = [ 9.5000,  -0.7500, pi-atan2(0.75,9.5)]';
eta_d_end = [0 0 pi-atan2(0.75,9.5)]';

Eta_d=eta_d_0*(1-theta)+ theta*eta_d_end;

 angle_tan=Eta_d(3);   
 angle_new =k*2*pi+angle_tan;
 if ((Eta_d(3)-angle_new)^2>pi^2/4)
     disp(Eta_d(3))
     disp(angle_new)
     if (angle_new>Eta_d(3))
         k=k-1;
         angle_new=k*2*pi+angle_tan;
     elseif (angle_new<Eta_d(3))
        k=k+1;
        angle_new=k*2*pi+angle_tan;
     end
 end
 
Eta_d(3)=angle_new;

Eta_dtheta=-eta_d_0+eta_d_end;
Eta_dtheta2=[0 0 0]';
v_s=u_d/norm([9.5 0.75]);
v_stheta=0;
v_st=u_ddot/norm([9.5 0.75]);
theta_dot=v_s;
theta=euler2(theta_dot,theta,h); %% Integrating with euler
theta_r=theta;
end
