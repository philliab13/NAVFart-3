function [n_c]= open_loop_speed_control(U_ref)

    Dia = 3.3;              % propeller diameter (m)
    rho = 1025; 
    KT = 0.6367;   
    Xudot = -8.9830e5;
    T1 = 20; 

    m = 17.0677e6;   
    Xu = -(m-Xudot)/T1;

    T_d=U_ref*Xu/(0.05-1);
    n_c    = sqrt( max(T_d,0) / (rho*Dia^4*KT) );
    


    end