function matrix = D_c(V)

      global X_u_c Y_v_c Y_r_c N_v_c N_r_c;

      D_L=[-X_u_c    0    0;
            0     -Y_v_c -Y_r_c;
            0     -N_v_c -N_r_c];
      matrix=D_L+D_NL(V);

end

function matrix = D_NL(V)

       global X_uu_c Y_vv_c Y_vr_c Y_rv_c Y_rr_c N_vv_c N_vr_c N_rv_c N_rr_c X_uuu_c


    u=V(1);
    v=V(2);
    r=V(3);
     matrix=[-X_uu_c*abs(u)-X_uuu_c*u^2 0 0;
            0 -Y_vv_c*abs(v)-Y_rv_c*abs(r) -Y_vr_c*abs(v)-Y_rr_c*abs(r);
            0 -N_vv_c*abs(v)-N_rv_c*abs(r) -N_vr_c*abs(v)-N_rr_c*abs(r)];
end