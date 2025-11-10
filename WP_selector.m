function [x_t, y_t, x_ref, y_ref, last_wp] = WP_selector(x, y, varargin)
% WP_SELECTOR
%  x_ref,y_ref : current/reference waypoint (k)
%  x_t,y_t     : target waypoint (k+1)  [or hold at last if none]
%  last_wp     : 1 iff you're inside the acceptance radius of the final WP

persistent WP k

if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1},'reset')
    WP = []; k = [];
end

R_accept = 2*161;  % acceptance radius (meters)

% Load once
if isempty(WP)
    S  = load('WP.mat','WP'); 
    WP = S.WP;                 % expect 2×N (rows: x;y)
    if size(WP,1) ~= 2 && size(WP,2) == 2
        WP = WP.';             % fix if provided as N×2
    end
    k = 1;                     % start at first waypoint
end

N = size(WP,2);
k = min(max(k,1), N);          % clamp for safety

% Distance^2 to current reference (k)
dx = WP(1,k) - x; 
dy = WP(2,k) - y; 
if (k < N) && (dx*dx + dy*dy <= R_accept^2)
    k = k + 1;                 % advance when inside acceptance of (k)
end

% Outputs: reference (k) and target (k+1 or hold)
x_ref = WP(1,k); 
y_ref = WP(2,k);
if k < N
    x_t = WP(1,k+1); 
    y_t = WP(2,k+1);
else
    x_t = x_ref; 
    y_t = y_ref;              % hold last WP as target
end

% Final waypoint flag: inside acceptance of the last WP
dxF = WP(1,N) - x; 
dyF = WP(2,N) - y;
last_wp = double(k == N && (dxF*dxF + dyF*dyF <= R_accept^2));
end
