function [phi,Omega_n] = Adaption_law_cl(Eta,nu,nu_dot,phi,Omega_n,tau,h)

persistent PHI OMEGA_N 

if (isempty(PHI)||isempty(OMEGA_N))
    PHI=0;
    OMEGA_N=0;
end

global M_c M_c_inv Gamma_phi Gamma_omega

R=Rot(Eta);

Y=M_c*nu_dot+C_c(nu)*nu+G_d(nu)-tau;
V=-Phi(nu)*phi+R'*Omega_n;
epsilon=Y-V;


phi_dot=@(phi_v,v)(-Gamma_phi*Phi(v)'*epsilon-PHI);
Omega_n_dot=@(omega_v,v)(Gamma_omega*R*epsilon-OMEGA_N);


PHI=PHI;%+Gamma_phi*Phi(nu)'*epsilon*0;
OMEGA_N=OMEGA_N-0.1*h*Gamma_omega*R'*epsilon;

phi=rk4(phi_dot,phi,nu,h);
Omega_n=rk4(Omega_n_dot,Omega_n,nu,h);
%Omega_n=Omega_n*0;



end

