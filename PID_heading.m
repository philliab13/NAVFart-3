function [delta_c]= PID_heading(e_psi, e_r, e_int)
        %Kp=195.6;
        %Kd=4073;
        %Ki=1.82;
        w_n=0.093;
        T=-99.47;
        K=-0.0049;

        Kp=w_n^2*T/K;
        Kd=(2*w_n*T-1)/K;
        Ki=Kp*0.1;

        
        delta_c=-(Kp*ssa(e_psi)+Kd*e_r+Ki*ssa(e_int));


    end