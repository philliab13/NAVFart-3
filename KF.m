function [x_pst, P_pst, x_prd, P_prd] = KF(x_prd, P_prd, Ad, Bd, Ed, Cd, Qd, Rd, psi_meas, rudder)
n = size(Ad, 1);

x_prd = Ad * x_prd + Bd * rudder;
P_prd = Ad * P_prd * Ad' + Ed * Qd * Ed';

K = P_prd * Cd' / (Cd * P_prd * Cd' + Rd);

y_tilde = psi_meas - Cd * x_prd;

x_pst = x_prd + K * y_tilde;
P_pst = (eye(n) - K * Cd) * P_prd * (eye(n) - K * Cd)' + K * Rd * K';

end






%{
function [x_hat,P_hat,x_prd,P_prd] = KF(x_prd,P_prd,Ad,Bd,Ed,Cd,Qd,Rd,psi_meas,rudder)
N = 1;

for i = 1 : N
    K = P_prd * Cd' * inv(Cd * P_prd * Cd' + Rd);
    IKC = eye(n) - K * Cd;

    y = psi_meas;

    x_hat = x_prd + K * (y - Cd * x_prd);
    P_hat = IKC * P_prd * IKC' + K * Rd * K';

    x_prd = Ad * x_hat + Bd * rudder;
    P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';
end
%}