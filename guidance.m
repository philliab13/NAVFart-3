function [chi_d] = guidance(e_yp,pi_p)
%LOS-Guidance - Takes in crosstrack and computes the chi_d

delta=600;
Kp=1/delta;


chi_d = pi_p-atan(Kp*e_yp);

end