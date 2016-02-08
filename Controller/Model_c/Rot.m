function matrix=Rot(N)
    yaw=N(3);
    matrix=[cos(yaw) -sin(yaw) 0;
            sin(yaw)  cos(yaw) 0;
            0         0        1]; 
end
