%test C(v)
CC=[];
model_parameters_c
for i=1:100
    
    v=[i i+1 i+2]';
    CC=[CC C_c(v)*v D_c(v)*v Rot(v)*v ];
    
end

disp('test ok!')