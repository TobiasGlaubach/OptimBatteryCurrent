function theta = penalty_function(u, sigma, method, p)
% penalty_function  implements different forms of objective functions.
%   C = penalty_function(u, sigma, method, p)
%       defaults: 
%           method = 'logbarrier'
%           p = 2
%   method = 'norm' for norm with parameter p
%   method = 'deadzone for Deadzone-linear penalty function
%   method = 'logbarrier' for Log barrier penalty function
%
%   p is the norm number and is only used if method=0

    if nargin <= 2
        method = 'logbarrier';
    end

    if nargin <= 3
        p = 2;
    end

    switch method

    case 'norm'
        theta = norm(u, p);
    case 'deadzone'
        theta = deadzone(u, sigma);
    case 'logbarrier'
        theta = log_barrier(u, sigma);
    otherwise
        error('method type not understood')
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
        theta = -sigma^2 * log(1-(u/sigma).^2);
    else
        theta = inf;
    end
end