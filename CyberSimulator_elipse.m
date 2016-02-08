    clear all
   % close all
    clc
    
    scenario='test';
    
    
    addpath('model/');
    addpath('Controller/');
     addpath('Controller/Model_c');
    addpath('gnc_mfiles/');
    %% model
    model_parameters
    %% Control parameters
    model_parameters_c
    
     %% Simulation options
    t_end = 500;
    h = 0.01;   %% Find good steplength

    time = 0:h:t_end;
    n=length(time);

    %% Initial conditions
    Eta_0 = [ 11.5000,  -0.7500*0, 1.5*pi]';
    nu_0  = [0.15, 0,0]';
    
    %% controler starter
    theta=0;
    u_d=0.15;
    u_ddot=0;
    phi=0*[-0.805 -7.250 -0.845 -3.45 0.130 -1.9 0.080 -0.750]';
    Omega_n_s=zeros(6,1);
    
     Eta_t0=[ 9.5000   -0.7500  pi]';
    Energy = 0;

    IAE_pos = 0;
    ISE_pos = 0;
    ITAE_pos = 0;
    IAEW_pos = 0;
    
    
    %% Controller input 
    
    %% Disturbance model

    
  steptime=10;
%   w=[3.5 3.5 0]';
%   Omega_r_s=kron([zeros(1,floor(steptime*n/t_end)) ones(1,n-floor(steptime*n/t_end))],w); %Prosess noice

    mu=0;
    sigma=0.3; 
    w=[sin(time/5).^2+cos(time/5); sin(time/5).^2+cos(time/5); zeros(1,n)];
    rand=normrnd(mu*ones(3,n),sigma*ones(3,n));
    Omega_r_s=w+rand;
   
   %Prosess noice
    
    
%% Initialize logger
log = Logger();

log.init(time, h);
log.add('Eta', 3);
log.add('Eta_t',3);
log.add('Eta_tilde',3);
log.add('Chi_tilde',1);
log.add('Nu', 3);
log.add('tau', 3);
log.add('tau_L', 3);
log.add('tau_U', 3);
log.add('phi', 8);
log.add('omega_n', 3);
log.add('Omega_n_s', 3);

log.add('u', 3);
log.add('dNu', 3);
log.add('Power', 1);
log.add('Energy', 1);
log.add('e', 1);
%log.add('IAE_pos', 1);
log.add('ISE_pos', 1);
log.add('ITAE_pos', 1);
log.add('IAEW_pos', 1);
log.add('IAEWWT',1); %integral squear error work ware and tear
log.add('IAEW',1); 
log.add('AE_pos', 1);
log.add('SE_pos', 1);
log.add('theta',1);

%% Simulatiom starter
    
Eta = Eta_0;
nu  = nu_0;
Eta_hat = Eta;
nu_hat = nu;
nu_dot=zeros(3,1);
Omega_n=zeros(3,1)  ;


tau = zeros(3,1);  
tau_old = zeros(3,1);  
ISE=0;
IAEWWT=0;
ABSE=0;
E=0;
WT=0;

Eta_t=Eta_t0;
for k=1:length(time)
    t=time(k);
    disp((t_end-t)/t_end) 
   %% controller Update
    [Eta_t,Eta_ttheta,Eta_ttheta2,v_s,v_stheta,v_st,theta]=Guidance_elip_fast(Eta_t, u_d, u_ddot,h);
    Eta_tdot=Eta_ttheta*v_s;
    Eta_tdotdot=Eta_ttheta2*v_s*v_s+Eta_ttheta*v_stheta*v_s;
    [ tau,phi,Omega_n,z_1,tau_L,tau_U,Chi_tilde ] = Adaptiv_controller_IandI(Eta,nu,nu_dot,phi,Omega_n,Eta_t,Eta_tdot,Eta_tdotdot,tau,h);
   %% omega
    Omega_n_s_dot=@(Omega_n_s,Omega_r_s)(A_d*Omega_n_s+B_d*Omega_r_s);
    Omega_n_s=rk4(Omega_n_s_dot,Omega_n_s,Omega_r_s(:,k),h); 
    
    %% Simulation dynamics
    Eta_dot=@(Eta_v,nu_v)(Rot(Eta_v)*nu_v);
    Eta_next=rk4(Eta_dot,Eta,nu,h);
    
    f=@(v,tu)(M_s\(Rot(Eta)'*Omega_n_s(1:3)+tu-C_s(v)*v-D_s(v)*v));   %%vi må se på hvordan få forstyrrelse inn i ligningen
    nu_dot=feval(f,nu,tau);
    %[phi,Omega_n] = Adaption_law_cl(Eta,nu,nu_dot,phi,Omega_n,tau,h); %%% CL
   % 
    
   % f=M_s\(Rot(Eta)'*W(:,k)+tau-C_s(nu)*nu-D_s(nu)*nu);   %%vi må se på hvordan få forstyrrelse inn i ligningen
   % nu_next=euler2(f,nu,h);
    nu_next=rk4(f,nu,tau,h); 
    %% Time update
    Eta=Eta_next;
    nu=nu_next;
    
    %% Preformance Metrics
    e=norm(Eta([1 2])-Eta_t([1 2]));
    e=(Eta([1 2])-Eta_t([1 2]))'*(Eta_t([2 1])/norm(Eta_t([2 1]))); %crostrack error
    
    trackpath_ort=cross([cos(Eta_t(3)) sin(Eta_t(3)) 0]',[0 0 -1]');
    e=(Eta([1 2])-Eta_t([1 2]))'*trackpath_ort(1:2);
    
    sqe=e^2;
    abse=abs(e);
    ABSE=euler2(abse,ABSE,h);
    P=abs(nu'*tau);
    E=euler2(P,E,h);
    wt=norm((tau-tau_old)/h);
    WT=euler2(wt,WT,h);

    ISE=euler2(e^2,ISE,h);
    IAEWWT=ABSE*E*WT;
    IAEW=ABSE*E;
    
    tau_old=tau;
  
    %% storing    
    log.store('Eta', Eta, k);
    log.store('Eta_t',Eta_t,k);
    log.store('Eta_tilde',z_1,k);
    log.store('Chi_tilde',Chi_tilde,k);
    
    log.store('Nu', nu, k);
    log.store('tau', tau,k);
    log.store('tau_L', tau_L,k);
    log.store('tau_U', tau_U,k);
    log.store('phi', phi,k);
    
    log.store('omega_n',Omega_n, k);
    log.store('Omega_n_s',Omega_n_s(1:3), k);
    
    log.store('ISE_pos',ISE,k);
    log.store('theta',theta,k);
    log.store('e',e,k);
    log.store('IAEWWT',IAEWWT,k);
    log.store('IAEW',IAEWWT,k);
        
end
%% PRINTS
% close all
%line='-.g';
line=':k';
%line='-.m';

figure(1)
North=log.get('Eta',1);
East=log.get('Eta',2);
Yaw=log.get('Eta',3);
North_d=log.get('Eta_t',1);
East_d=log.get('Eta_t',2);
Yaw_d=log.get('Eta_t',3);
hold on
plot(East,North,line)
plot(East_d,North_d,'r')
scatter(East(1),North(1),'xr')
%quiver(East(1:1500:n/1.5),North(1:1500:n/1.5),sin(Yaw(1:1500:n/1.5)),cos(Yaw(1:1500:n/1.5)),'b')
%quiver(East_d(1:1000:n),North_d(1:1000:n),sin(Yaw_d(1:1000:n)),cos(Yaw_d(1:1000:n)))
axis([-10 10 -5 15])
xlabel('east')
ylabel('north')
% mkdir('figures',scenario)
% saveas(gcf,['figures/' scenario '/trajectory'], 'png')


figure(2)
subplot(2,1,1);
hold on
plot(time, log.get('Chi_tilde', 1),line);
title('Chi_Tilde')

subplot(2,1,2);
hold on
plot(time, log.get('e', 1),line);
title('e')

% mkdir('figures',scenario)
% saveas(gcf,['figures/' scenario '/errorplot'], 'png')

figure(3)
koff=['Y|r|v';'Yr   ';'Y|v|r';'Y|r|r';'N|r|v';'Nr   ';'N|v|r';'N|r|r'];
koef_s=[-0.805 -7.250 -0.845 -3.45 0.130 -1.9 0.080 -0.750];
PPHI=log.get('phi');
for i=1:8
subplot(4,2,i)
hold on
plot(time,PPHI(i,:),line)
plot(time,ones(1,length(time))*koef_s(i),'r')
title(koff(i,:))
end
% mkdir('figures',scenario)
% saveas(gcf,['figures/' scenario '/hydro_estimates'], 'png')

figure(4) 
subplot(3,1,1)
hold on
plot(time, log.get('omega_n', 1),line)
plot(time, log.get('Omega_n_s', 1),'r')
legend('omega_n','Omega_n_s')
title('\omega estimate north')
subplot(3,1,2)
hold on
plot(time, log.get('omega_n', 2),line)
plot(time, log.get('Omega_n_s', 2),'r')
legend('omega_n','Omega_n_s')
title('\omega estimate east')
subplot(3,1,3)
hold on
plot(time, log.get('omega_n', 3),line)
plot(time, log.get('Omega_n_s', 3),'r')
legend('omega_n','Omega_n_s')
title('\omega estimate rot')
% mkdir('figures',scenario)
% saveas(gcf,['figures/' scenario '/disturbance_estimates'], 'png')

figure(5)
subplot(2,1,1)
hold on
plot(time, log.get('ISE_pos', 1),line)
title('ISE_p_o_s')

subplot(2,1,2)
hold on
plot(time, log.get('IAEWWT', 1),line)
title('IAEWWT')


figure(6)
subplot(3,1,1)
hold on
tau_L=sqrt(sum(log.get('tau_L').^2));
%TAU_L=[TAU_L;tau_L];
plot(time(500:n),tau_L(500:n),line)
subplot(3,1,2)

hold on
tau_U=sqrt(sum(log.get('tau_U').^2));
%TAU_U=[TAU_U;tau_U];
plot(time(500:n),tau_U(500:n),line)
subplot(3,1,3)
hold on
%tau=sqrt(sum((log.get('tau_L')+log.get('tau_U')).^2));
tau=sqrt(sum(log.get('tau').^2));
%tau=log.get('tau');
%tau=tau(3,:);
%TAU=[TAU;tau];
plot(time(500:n),tau(500:n),line)

% mkdir('figures',scenario)
% saveas(gcf,['figures/' scenario '/preformance_metrics'], 'png')
%% Prints
% close all
%%

% %%
% figure(1)
% North=log.get('Eta',1);
% East=log.get('Eta',2);
% Yaw=log.get('Eta',3);
% North_d=log.get('Eta_t',1);
% East_d=log.get('Eta_t',2);
% Yaw_d=log.get('Eta_t',3);
% subplot(2,1,1)
% hold on
% plot(East,North,'-')
% plot(East_d,North_d,'r')
% scatter(East(1),North(1),'xr')
% %quiver(East(1:1000:n),North(1:1000:n),sin(Yaw(1:1000:n))*10,cos(Yaw(1:1000:n))*10)
% %quiver(East_d(1:1000:n),North_d(1:1000:n),sin(Yaw_d(1:1000:n))*10,cos(Yaw_d(1:1000:n))*10)
% axis([-50 50 -20 20])
% xlabel('east')
% ylabel('north')
% subplot(2,1,2)
% hold on
% plot(time,Yaw,'*-')
% plot(time,log.get('Eta_t',3),'r')
% 
% figure(2); clf
% 
% subplot(5,1,1);
% title('x');
% hold on
% plot(time, log.get('Eta', 1));
% plot(time, log.get('Eta_t',1),'r');
% hold off
% 
% subplot(5,1,2);
% title('y');
% hold on
% 
% plot(time, log.get('Eta', 2));
% plot(time, log.get('Eta_t',2),'r');
% hold off
% 
% subplot(5,1,3);
% title('psi');
% hold on
% plot(time, log.get('Eta', 3));
% plot(time, log.get('Eta_t',3),'r');
% hold off
% 
% subplot(5,1,4);
% hold on
% plot(time, log.get('Chi_tilde', 1));
% title('Chi_Tilde')
% 
% subplot(5,1,5);
% hold on
% plot(time, log.get('e', 1));
% title('e')
% 
% figure(3); clf
% 
% subplot(3,1,1);
% title('u');
% hold on
% plot(time, log.get('Nu', 1));
% hold off
% 
% subplot(3,1,2);
% title('v');
% hold on
% plot(time, log.get('Nu', 2));
% hold off
% 
% subplot(3,1,3);
% title('r');
% hold on
% plot(time, log.get('Nu', 3));
% hold off
% figure(4) 
% subplot(3,2,1)
% plot(time, log.get('tau', 1));
% title('tau_u')
% legend('tau_u')
% subplot(3,2,3)
% plot(time, log.get('tau', 2));
% title('tau_v')
% %legend('tau_u')
% subplot(3,2,5)
% plot(time, log.get('tau', 3));
% title('tau_r')
% subplot(3,2,2)
% plot(time, log.get('Eta_tilde', 1));
% title('Eta_tilde_x')
% subplot(3,2,4)
% plot(time, log.get('Eta_tilde', 2));
% title('Eta_tilde_y')
% subplot(3,2,6)
% plot(time, log.get('Eta_tilde', 3));
% title('Eta_tilde_yaw')
% 
% figure(5) 
% subplot(3,2,1)
% plot(time, log.get('tau_L', 1));
% title('tau_uL')
% %legend('tau_uL')
% subplot(3,2,3)
% plot(time, log.get('tau_L', 2));
% title('tau_vL')
% %legend('tau_vL')
% subplot(3,2,5)
% plot(time, log.get('tau_L', 3));
% title('tau_rL')
% %legend('tau_rL')
% subplot(3,2,2)
% plot(time, log.get('tau_U', 1));
% title('tau_uU')
% 
% %legend('tau_uU')
% subplot(3,2,4)
% plot(time, log.get('tau_U', 2));
% title('tau_vU')
% %legend('tau_vU')
% subplot(3,2,6)
% plot(time, log.get('tau_U', 3));
% title('tau_rU')
% %legend('tau_rU')
% 
% figure(6)
% koff=['Y|r|v';'Yr   ';'Y|v|r';'Y|r|r';'N|r|v';'Nr   ';'N|v|r';'N|r|r'];
% koef_s=[-0.805 -7.250 -0.845 -3.45 0.130 -1.9 0.080 -0.750];
% PPHI=log.get('phi');
% for i=1:8
% subplot(4,2,i)
% hold on
% plot(time,PPHI(i,:))
% plot(time,ones(1,length(time))*koef_s(i),'r')
% title(koff(i,:))
% end
% 
% figure(7)
% subplot(2,1,1)
% plot(time, log.get('ISE_pos', 1))
% subplot(2,1,2)
% plot(time, log.get('ISEWWT', 1))
save(['LOGS/' scenario],'log')
