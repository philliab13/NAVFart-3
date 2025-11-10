function [chi_d] = ILOS_guidance(pi_p,e_yp,e_yint)
%UNTITLED Summary of this function goes here


delta=600;
Kp=1/delta;
Kappa=0.7;
Ki=Kappa*Kp;


chi_d = pi_p-atan(Kp*e_yp+Ki*e_yint);
end