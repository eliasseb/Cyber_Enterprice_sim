%Boddy class, contains the functions and variables for the boddy cinematics
%

classdef Boddy
   properties (GetAccess = 'public', SetAccess = 'private')
      u=0;
      v=0;
      r=0;
   end
   methods
      function obj = Boddy(initvalue)
           model_parameters
           obj.u=initvalue(1);
           obj.v=initvalue(2);
           obj.r=initvalue(3);
      end
      function state = V(obj) %returnerer "farts vektoren"
           state=[obj.u obj.v obj.r]';
      end
      function matrix=C_RB(obj)
        model_parameters
        matrix=[0                 0         -m*x_g*obj.r+obj.v;
               0                 0          m*obj.u;
               m*x_g*obj.r+obj.v m*obj.u    0]; 
      end
      function matrix=C_A(obj)
        model_parameters
        matrix=[0                                      0                Y_vdot*obj.v+0.5*(N_vdot+Y_rdot)*obj.r;
               0                                       0                X_udot*obj.u;
               -Y_vdot*obj.v+0.5*(N_vdot+Y_rdot)*obj.r -X_udot*obj.u    0]; 
      end      
      function matrix=C(obj) %%Coreolis matrix
        model_parameters
        matrix=C_RB(obj)+C_A(obj);
      end
   end
end
