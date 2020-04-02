function theta = penalty_function(u, sigma, method, p)
% penalty_function  implements different forms of objective functions.
%   C = penalty_function(u, sigma, method, p)
%       defaults: 
%           method = 2
%           p = 2
%   method = 0 for norm with parameter p
%   method = 1 for Deadzone-linear penalty function
%   method = 2 for Log barrier penalty function
%
%   p is the norm number and is only used if method=0

    if nargin <= 2
        method = 2;
    end

    if nargin <= 3
        p = 2;
    end

    switch method

    case 0
        theta = norm(u, p);
    case 1
        theta = deadzone(u, sigma);
    case 2
        theta = log_barrier(u, sigma);
    end

    return
end

function theta = deadzone(u, sigma)

    if abs(u) < sigma
        theta = 0;
    else
        theta = abs(u);
    end
end

function theta = log_barrier(u, sigma)
    if abs(u) < sigma
        theta = -sigma^2 * log10(1-(u/sigma)^2);
    else
        theta = inf;
    end
end