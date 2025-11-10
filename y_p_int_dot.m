function [dot] = y_p_int_dot(e_y_p,e_p_int)
%UNTITLED2 Summary of this function goes here

delta=600;
kappa=0.7;

dot=(delta*e_y_p)/(delta^2+(e_y_p+kappa*e_p_int)^2);

end