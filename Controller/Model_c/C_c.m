%coreolis matrix
%
    function matrix=C_c(V) %%Coreolis matrix
        matrix=C_RB_c(V)+C_A_c(V);
    end
      function matrix=C_RB_c(V)
            global x_g_c m_c
            u=V(1);
            v=V(2);
            r=V(3);
            matrix=[0        0       -m_c*(x_g_c*r+v);
                    0        0        m_c*u;
                    m_c*(x_g_c*r+v) -m_c*u    0]; 
      end
      function matrix=C_A_c(V)
            global X_udot_c Y_vdot_c Y_rdot_c N_vdot_c

            u=V(1);
            v=V(2);
            r=V(3);
        matrix=[0                                      0              Y_vdot_c*v+0.5*(N_vdot_c+Y_rdot_c)*r;
               0                                       0             -X_udot_c*u;
               -Y_vdot_c*v-0.5*(N_vdot_c+Y_rdot_c)*r         X_udot_c*u       0]; 
      end      

     
   

