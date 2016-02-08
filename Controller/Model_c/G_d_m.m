function [ d ] = G_d_m(nu)

global X_u_c Y_v_c N_v_c ;
    v=[nu(1) nu(2) 0]';

      D_L_uv=[-X_u_c    0  0 ;
            0     -Y_v_c 0;
            0     -N_v_c 0];
      d=(D_L_uv+D_NL_uv(v));

end

function matrix = D_NL_uv(V)

       global X_uu_c Y_vv_c N_vv_c X_uuu_c


    u=V(1);
    v=V(2);
    
     matrix=[-X_uu_c*abs(u)-X_uuu_c*u^2 0 0;
            0 -Y_vv_c*abs(v) 0;
            0 -N_vv_c*abs(v) 0];
end
