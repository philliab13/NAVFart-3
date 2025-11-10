   function [xd_dot] =ref_model(xd,psi_ref)
        wn=0.03;
        zeta=1;
        
        A_d=[0 1 0;
             0 0 1; 
             -wn^3 -(2*zeta+1)*wn^2 -(2*zeta+1)*wn;];
        B_d=[0;
           0;
           wn^3;];

        


        xd_dot=A_d*xd +B_d*psi_ref;
    end