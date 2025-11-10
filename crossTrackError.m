function [e_y, pi_p] = crossTrackError(xkp1, ykp1, xk, yk, x, y)
dx = xkp1 - xk;  dy = ykp1 - yk;
L2 = dx*dx + dy*dy;
if L2 < 1e-6
    % Degenerate segment: treat as point target
    pi_p = atan2(ykp1 - y, xkp1 - x);
    e_y  = 0;
    return
end
pi_p = atan2(dy, dx);
% cross-track error in path frame
e_y  = -sin(pi_p) * (x - xk) + cos(pi_p) * (y - yk);
end