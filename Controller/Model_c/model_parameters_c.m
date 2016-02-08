%% Av elias Bjørne, for prosjektoppgave adaptiv ulinær regulator
%% Modellen er av Cybership II, basert på parametere fra Roger sjetne 2004
%% Bruksmåte 


%% Parametre for modellen 
   global  g x_g_c m_c I_z_c X_udot_c Y_vdot_c Y_rdot_c N_vdot_c N_rdot_c X_u_c Y_v_c Y_r_c N_v_c N_r_c X_uu_c Y_vv_c Y_vr_c Y_rv_c Y_rr_c N_vv_c N_vr_c N_rv_c N_rr_c X_uuu_c M_c M_RB_c M_A_c D_L_c M_c_inv A_d B_d
%% Parameters for controler
    global K_p K_d Gamma_phi Gamma_omega
    
   K_p=diag([0.2,1,0.1]);
   K_d=diag([5,7,6]);
%    K_p=diag([0.4,2,0.2]);
%    K_d=diag([5,7,6]);
%    

%    K_p=diag([4.3996,2.5876,2.7812]);
%    K_d=diag([12.0990,19.9655,16.7353]);

%    K_p=diag([0.9780,4.9829 ,2.6622]);
%    K_d=diag([4.6093,17.4250,8.7047]);
   
   Gamma_phi=diag([8 4 8 8 8 4 8 8]*10);
   Gamma_omega=diag([1 1 0.1]);

   %Generic parameters
    g=9.81;
    %% Controll model parameters
    m_c = 23.8;
    I_z_c = 1.760;
    x_g_c = 0.046;

    X_udot_c = -2.0;
    Y_vdot_c = -10.0;
    Y_rdot_c = -0.0;
    N_vdot_c = -0.0;
    N_rdot_c = -1.0;

    X_u_c = -0.72253;
    Y_v_c = -0.88965;
    Y_r_c = -7.250;
    N_v_c = 0.03130;
    N_r_c = -1.900;

    X_uu_c = -1.32742;
    Y_vv_c = -36.47287;
    Y_vr_c = -0.845;
    Y_rv_c = -0.805;
    Y_rr_c = -3.450;
    N_vv_c = 3.95645;
    N_vr_c = 0.080;
    N_rv_c = 0.130;
    N_rr_c = -0.750;

    X_uuu_c = -5.86643;
%% Matrices
  
    
    M_RB_c=[m_c     0      0;
          0     m_c      m_c*x_g_c;
          0     m_c*x_g_c  I_z_c];
    M_A_c=[-X_udot_c    0          0;
            0          -Y_vdot_c  -Y_rdot_c;
            0          -N_vdot_c  -N_rdot_c];     
    M_c=M_RB_c+M_A_c;
    
    D_L_c=[-X_u_c    0    0;
         0      -Y_v_c -Y_r_c;
         0      -N_v_c -N_r_c];
     
     M_c_inv = [   0.0388         0         0      ;
                0              0.0300   -0.0119 ;
                0             -0.0119    0.3670 ];
            
            
    %% Refrence Model
    OMEGA=diag([1 1 1]);
    LAMBDA=diag([1 1 2]);
    
    A_d=[ zeros(3,3) eye(3)
         -OMEGA^2 -2*LAMBDA*OMEGA];
     
     B_d=[zeros(3,3);OMEGA^2];
     
     