function [x_t, y_t, x_ref, y_ref, last_wp] = WP_selector(xN, yE, varargin)
% xN,yE in North-East frame

persistent WP k
if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1},'reset')
    WP = []; k = [];
end

R_accept = 2*161; % m

% Load once
if isempty(WP)
    S = load('WP.mat','WP');
    WP = S.WP;                     % expect 2×N (rows: [x;y]=[North;East])
    if size(WP,1) ~= 2 && size(WP,2) == 2
        WP = WP.';                 % if given as N×2, transpose
    end
    k = 1;
end

N = size(WP,2);
k = min(max(k,1), N);

pos = [xN; yE];

% Advance while conditions are met (handles big steps or very large R)
while k < N
    wk  = WP(:,k);
    wk1 = WP(:,k+1);
    v   = wk1 - wk;
    L2  = v.'*v;

    % Degenerate segment: skip it
    if L2 < 1e-8
        k = k + 1;
        continue
    end

    % Tests
    near_curr = sum((pos - wk ).^2) <= R_accept^2;
    near_next = sum((pos - wk1).^2) <= R_accept^2;

    % Half-plane at w_{k+1} (switch when you've "passed" k+1 along v)
    passed_next = (pos - wk1).' * v >= 0;

    if near_curr || near_next || passed_next
        k = k + 1;         % move to next reference
    else
        break
    end
end

% Outputs
x_ref = WP(1,k); y_ref = WP(2,k);
if k < N
    x_t = WP(1,k+1); y_t = WP(2,k+1);
else
    x_t = x_ref;   y_t = y_ref;    % hold last
end

% Final waypoint flag
last_wp = double( k==N && sum((pos - WP(:,N)).^2) <= R_accept^2 );
end
