function matrix = D_s(V)
    X_u = -0.655;
    Y_v = -1.330;
    Y_r = -0;
    N_v = 0;
    N_r = -1.900;
    
    D_L=[-X_u    0    0;
          0     -Y_v -Y_r;
          0     -N_v -N_r];
    matrix=D_L+D_NL(V);

end

function matrix = D_NL(V)
    X_uu = 0.355;
    Y_vv = -2.776;
    Y_vr = -0;
    Y_rv = -0;
    Y_rr = -0.845;
    N_vv = 0.805;
    N_vr = 0;
    N_rv = 0;
    N_rr = -0.750;
    X_uuu = -0;

    u=V(1);
    v=V(2);
    r=V(3);
     matrix=[-X_uu*abs(u)-X_uuu*u^2 0 0;
            0 -Y_vv*abs(v)-Y_rv*abs(r) -Y_vr*abs(v)-Y_rr*abs(r);
            0 -N_vv*abs(v)-N_rv*abs(r) -N_vr*abs(v)-N_rr*abs(r)];
end