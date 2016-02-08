function [ tau,phi,Omega_n,z_1,tau_L,tau_U,Chi_tilde ] = Adaptiv_controller_cl(Eta,nu,nu_dot,phi,Omega_n,Eta_d,Eta_ddot,Eta_ddotdot,tau,h)
%  [ tau,theta ] = Adaptiv_controller(Eta,nu,nu_dot,phi,theta_t,Eta_d,Eta_dtheta,Eta_dtheta2,v_s,v_stheta,v_st,h)
%UNTITLED6 Summary of this function goes here Skjetne
%   Detailed explanation goes here
%kan være den må løses på en smartere måte
%include Model_c to path

 global K_d Gamma_phi Gamma_omega M_c M_c_inv K_p
 
R=Rot(Eta);
  
%% controlsignaler
Eta_tilde=Eta-Eta_d;
Chi_tilde=Eta(3)-Eta_d(3);
Chi_tilde=mod(pi/2+Chi_tilde,pi)-pi/2; %%%%%%%%%%%%%%%%%%%%%%

%% New Controller
Eta_tilde_dot=R*nu-Eta_ddot;

z_1=R'*(Eta_tilde);
z_1_dot=S(nu)'*z_1+Rot_dot(nu)'*Eta_tilde + R'*Eta_tilde_dot;

Alpha_1=(-K_p*z_1+(R')*Eta_ddot);
Alpha_1dot=-K_p*z_1_dot+Rot_dot(nu)'*Eta_ddot*h + R'*Eta_ddotdot*h; %%%%%%%%%%% Tatt bort 

z_2=nu-Alpha_1;

%% estimering av Phi

phi=phi;
%% estimering av omega_n
Omega_n=Omega_n;

%% (nu+nu_dot*h)
% TAU 
%%%%%%%%%%%%%%tau=C_c(nu)*Alpha_1+G_d_m(nu)*Alpha_1+Phi(nu)*phi+M_c*Alpha_1dot-z_1-K_d*z_2-R'*omega_n;
%tau=C_c(nu)*nu+M_c*Alpha_1dot-z_1-K_d*z_2-R'*omega_n;
tau=C_c(nu)*nu+G_d(nu)+Phi(nu)*phi+M_c*Alpha_1dot-z_1-K_d*z_2-R'*Omega_n;


%%
tau_L=-z_1-K_d*z_2+M_c*Alpha_1dot;
tau_U=C_c(nu)*nu+G_d(nu)+Phi(nu)*phi-R'*Omega_n;


end

