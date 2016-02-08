%% Av elias Bjørne, for prosjektoppgave adaptiv ulinær regulator
%% Modellen er av Cybership II, basert på parametere fra Roger sjetne 2004
%% Bruksmåte 


%% Parametre for modellen 
   global  g x_g m I_z X_udot Y_vdot Y_rdot N_vdot N_rdot X_u Y_v Y_r N_v N_r X_uu Y_vv Y_vr Y_rv Y_rr N_vv N_vr N_rv N_rr X_uuu M_s M_RB M_A D_L
    %Generic parameters
    g=9.81;
    % Parameters ship model
    x_g = 0.0375;
    m = 14.790;
    I_z = 1.760;

    X_udot = -2.0;
    Y_vdot = -10.0;
    Y_rdot = -0.0;
    N_vdot = -0.0;
    N_rdot = -1.0;

    X_u = -0.655;
    Y_v = -1.330;
    Y_r = -0;
    N_v = 0;
    N_r = -1.900;

    X_uu = 0.355;
    Y_vv = -2.776;
    Y_vr = -0;
    Y_rv = -0;
    Y_rr = -0.845;
    N_vv = 0.805;
    N_vr = 0;
    N_rv = 0;
    N_rr = -0.750;
    
    X_uuu = -0; %?????????
%% Matrices
  
    
    M_RB=[m     0      0;
          0     m      m*x_g;
          0     m*x_g  I_z];
    M_A=[-X_udot    0        0;
         0          -Y_vdot  -Y_rdot;
         0          -N_vdot  -N_rdot];     
    M_s=M_RB+M_A;
    
    D_L=[-X_u    0    0;
         0      -Y_v -Y_r;
         0      -N_v -N_r];
    disp('param_laded')