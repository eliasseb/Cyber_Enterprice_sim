%Boddy class, contains the functions and variables for the boddy cinematics
%
    function matrix=C_s(V) %%Coreolis matrix
        matrix=C_RB(V)+C_A(V);
    end
      function matrix=C_RB(V)
       x_g = 0.0375;
        m = 14.790;
            
            u=V(1);
            v=V(2);
            r=V(3);
            matrix=[0        0       -m*(x_g*r+v);
                    0        0        m*u;
                    m*(x_g*r+v) -m*u    0]; 
      end
      function matrix=C_A(V)
        X_udot = -2.0;
        Y_vdot = -10.0;
        Y_rdot = -0.0;
        N_vdot = -0.0;


            u=V(1);
            v=V(2);
            r=V(3);
        matrix=[0                                      0              Y_vdot*v+0.5*(N_vdot+Y_rdot)*r;
               0                                       0             -X_udot*u;
               -Y_vdot*v-0.5*(N_vdot+Y_rdot)*r         X_udot*u       0]; 
      end      

     
   

